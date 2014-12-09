/****************************************************************
 unveil-eval-path.C
 An HMM genefinder based on Salzberg's VEIL model.  This version
 evaluates a given gene model.

 bmajoros@tigr.org
 Feb. 2003
 ****************************************************************/

#include <string>
#include <iostream.h>
#include "tigrstl/CommandLine.H"
#include "tigrstl/SetPlus.H"
#include "tigrstl/Regex.H"
#include "tigrstl/VectorPlus.H"
#include "tigrstl/MapPlus.H"
#include "tigrstl/GffReader.H"
#include "tigrstl/GffTranscript.H"
#include "tigrstl/GffExon.H"
#include "tigrstl/GffFeature.H"
#include <math.h>

#include "FastViterbi.H"
#include "FastaReader.H"

enum ExonTypes
  {
    INITIAL_EXON,
    INTERNAL_EXON,
    FINAL_EXON,
    SINGLE_EXON
  };
ostream &operator<<(ostream &,ExonTypes &);


class Exon
{
  int begin, end, frame, transcriptId;
  ExonTypes exonType;
  StringPlus axisId;
  char strand;
public:
  Exon(int begin,int end,int frame,int transcriptId,
       ExonTypes exonType,StringPlus axisId,char strand)
    : begin(begin), end(end), frame(frame), transcriptId(transcriptId),
      exonType(exonType), axisId(axisId), strand(strand) {}
  void printOn(ostream &);
  int getTranscriptId() const {return transcriptId;}
  ExonTypes getExonType() const {return exonType;}
  void setExonType(ExonTypes t) {exonType=t;}
  char getStrand() const {return strand;}
};
ostream &operator<<(ostream &,Exon &);


class Application
{
  SetPlus<int> startCodonStates, stopCodonStates, exonStates,
    upstreamIntergenicStates, downstreamIntergenicStates, frameshiftStates;
  Regex deflineRegex, donorRegex, acceptorRegex;
  vector<HMM*> submodels;
  int upstreamModelId, downstreamModelId, strandDelta, numSymbols;
  vector<int> submodelOrigins; // state ids in composite HMM
  int numSubmodels;
  VectorPlus<Exon*> exons;
  MapPlus<int,int> exonsPerGene;

  void copyIn(HMM &submodel,HMM &compositeHMM,int delta);
  void addToSet(SetPlus<int> &,int fromState,int numStates);
  void loadSubmodels(const StringPlus &filename);
  HMM *loadHMM(const StringPlus &);
  void linkModels(int fromModelId,int fromDelta,
		  int toModelId,int toDelta,
		  HMM &compositeHMM,double metaP);
  void adjustExonTypes();
  void produceOutput();
  double getPathScore(vector<int> &path,Sequence &,HMM &);
  void getCodingStates(SetPlus<int> &codingStates);
  BitMapSet *getCodingPositions(VectorPlus<GffTranscript*> &transcripts,
				 int seqLen);
public:
  Application();
  int main(int argc,char *argv[]);
};


int main(int argc,char *argv[])
  {
    try
      {
	Application app;
	return app.main(argc,argv);
      }
    catch(const char *p)
      {
	cerr << p << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



Application::Application()
  : deflineRegex(">\\s*(\\S+)"), 
    donorRegex("donor-(\\d+)"), acceptorRegex("acceptor-(\\d+)")
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    CommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=4)
      throw string("unveil-eval-path <*.hmm> <submodels.txt> <*.fasta> <*.gff>");
    StringPlus hmmFilename=cmd.arg(0);
    StringPlus submodelsFilename=cmd.arg(1);
    StringPlus fastaFilename=cmd.arg(2);
    StringPlus gffFilename=cmd.arg(3);

    loadSubmodels(submodelsFilename);
    HMM &hmm=*loadHMM(hmmFilename);
    Alphabet &alphabet=hmm.getAlphabet();
    FastViterbi viterbi(hmm);

    FastaReader fastaReader(fastaFilename);
    StringPlus defline, sequence;
    fastaReader.nextSequence(defline,sequence);
    if(!deflineRegex.search(defline))
      throw StringPlus("bad defline: ")+defline;
    StringPlus axisId=deflineRegex[1];
    Sequence seq(sequence,alphabet);
    int seqLen=seq.getLength();

    SetPlus<int> codingStates;
    getCodingStates(codingStates);

    // Load GFF
    GffReader gffReader(gffFilename);
    VectorPlus<GffTranscript*> *transcripts=gffReader.loadTranscripts();
    BitMapSet *codingPositions=getCodingPositions(*transcripts,seqLen);

    vector<int> &path=
      *viterbi.getPath_Masked(seq,*codingPositions,codingStates,strandDelta);

    double pathScore=getPathScore(path,seq,hmm);
    cout << "SCORE="<<pathScore << endl;

    return 0;
  }



void Application::getCodingStates(SetPlus<int> &codingStates)
{
  codingStates+=startCodonStates;
  codingStates+=stopCodonStates;
  codingStates+=exonStates;

  // ### HACK!!!!!!
  for(int i=114 ; i<=119 ; ++i) codingStates+=i;
  codingStates+=131;
  // ###
}



BitMapSet *Application::getCodingPositions(VectorPlus<GffTranscript*> 
					   &transcripts,
					   int seqLen)
{
  BitMapSet &codingPositions=*new BitMapSet(seqLen);

  int n=transcripts.size();
  for(int i=0 ; i<n ; ++i)
    {
      GffTranscript &transcript=*transcripts[i];
      int n=transcript.getNumExons();
      for(int i=0 ; i<n ; ++i)
	{
	  GffExon &exon=transcript.getIthExon(i);
	  int begin=exon.getBegin();
	  int end=exon.getEnd();
	  for(int i=begin ; i<end ; ++i)
	    codingPositions.addMember(i);
	}
    }

  return &codingPositions;
}



void Application::loadSubmodels(const StringPlus &filename)
{
  submodelOrigins.push_back(0);
  int nextFreeState=1;
  File file(filename);
  numSubmodels=0;
  while(!file.eof())
    {
      StringPlus line=file.readLine();
      if(file.eof()) break;
      vector<StringPlus> &fields=*line.getFields("= \t\n[]");
      int numFields=fields.size();
      if(numFields==1)
	throw StringPlus("Syntax error in submodel list file: ")+line;
      if(numFields>0)
	{
	  int modelId=fields[0].asInt();
	  StringPlus secondField=fields[1];
	  StringPlus modelFilename=fields[numFields-1];
	  HMM *submodel=new HMM(modelFilename);
	  while(modelId>=submodels.size()) submodels.push_back(NULL);
	  submodels[modelId]=submodel;
	  int submodelStates=submodel->getNumStates()-1;
	  ++numSubmodels;

	  if(secondField=="frameshift")
	    {
	      addToSet(frameshiftStates,nextFreeState,submodelStates);
	      addToSet(exonStates,nextFreeState,submodelStates);
	    }
	  else if(secondField=="start") 
	    {
	      addToSet(startCodonStates,nextFreeState,submodelStates);
	      addToSet(exonStates,nextFreeState,submodelStates);
	    }
	  else if(secondField=="stop")
	    {
	      addToSet(stopCodonStates,nextFreeState,submodelStates);
	      addToSet(exonStates,nextFreeState,submodelStates);
	    }
	  else if(secondField=="exon")
	    addToSet(exonStates,nextFreeState,submodelStates);
	  else if(secondField=="upstream-intergenic")
	    {
	      addToSet(upstreamIntergenicStates,nextFreeState,submodelStates);
	      upstreamModelId=modelId;
	    }
	  else if(secondField=="downstream-intergenic")
	    {
	      addToSet(downstreamIntergenicStates,nextFreeState,
		       submodelStates);
	      downstreamModelId=modelId;
	    }
	  else if(donorRegex.match(secondField))
	    {
	      // The donor model "steals" some number of exon nucleotides
	      int stolenDonorStates=donorRegex[1].asInt();
	      addToSet(exonStates,nextFreeState,stolenDonorStates);
	    }
	  else if(acceptorRegex.match(secondField))
	    {
	      // The acceptor model "steals" some number of exon nucleotides
	      int stolen=acceptorRegex[1].asInt();
	      addToSet(exonStates,nextFreeState+submodelStates-stolen,stolen);
	    }

	  submodelOrigins.push_back(nextFreeState);
	  nextFreeState+=submodelStates;
	}
      delete &fields;
    }
}



void Application::addToSet(SetPlus<int> &theSet,int fromState,int numStates)
{
  int end=fromState+numStates;
  for(int i=fromState ; i<end ; ++i)
    theSet.insert(i);
}



void Application::copyIn(HMM &submodel,HMM &compositeHMM,int delta)
{
  // Copy transitions
  int numStates=submodel.getNumStates();
  for(int i=1 ; i<numStates ; ++i)
    {
      for(int j=1 ; j<numStates ; ++j)
	compositeHMM.setTransProb(i+delta,j+delta,
				  submodel.getTransProb(i,j));
      compositeHMM.setTransProb(i+delta,0,
				submodel.getTransProb(i,0));
    }
  
  // Copy emissions
  for(int i=1 ; i<numStates ; ++i)
    for(int j=0 ; j<numSymbols ; ++j)
      compositeHMM.setEmitProb(i+delta,Symbol(j),
			       submodel.getEmitProb(i,Symbol(j)));
}



void Application::linkModels(int fromModelId,int fromStrandDelta,
			     int toModelId,int toStrandDelta,
			     HMM &compositeHMM,double metaP)
{
  /*
    This method takes every state in the "from" model having a
    transition to zero and redirects that transition so that
    it points to the state(s) that the zero in the "to" model
    can transition to.
   */

  HMM &fromModel=*submodels[fromModelId];
  HMM &toModel=*submodels[toModelId];
  int fromModelDelta=submodelOrigins[fromModelId]-1;
  int toModelDelta=submodelOrigins[toModelId]-1;
  int numToModelStates=toModel.getNumStates();
  int numFromModelStates=fromModel.getNumStates();
  int numAllStates=compositeHMM.getNumStates();

  for(int fromState=1 ; fromState<numFromModelStates ; ++fromState)
    {
      if(!fromStrandDelta && !fromModel.transExists(fromState,0)) continue;
      if(fromStrandDelta && !fromModel.transExists(0,fromState)) continue;
      for(int i=1 ; i<numAllStates ; ++i)
	{
	  if(i>=fromModelDelta+fromStrandDelta+1 && 
	     i<fromModelDelta+fromStrandDelta+numFromModelStates)
	    continue;
	  int adjustedFrom=fromState+fromModelDelta+fromStrandDelta;
	  double oldP=compositeHMM.getTransProb(adjustedFrom,i);
	  compositeHMM.setTransProb(adjustedFrom,i,oldP/2);
	}
      double fromP=fromModel.getTransProb(fromState,0);
      for(int toState=1 ; toState<numToModelStates ; ++toState)
	{
	  if(!toStrandDelta && !toModel.transExists(0,toState)) continue;
	  if(toStrandDelta && !toModel.transExists(toState,0)) continue;
	  double toP=toModel.getTransProb(0,toState);
	  double newP=fromP*metaP*toP;
	  if(compositeHMM.transExists(fromState+fromModelDelta+fromStrandDelta,
				      toState+toModelDelta+toStrandDelta))
	    throw StringPlus("ERROR: transition already exists!");
	  compositeHMM.setTransProb(fromState+fromModelDelta+fromStrandDelta,
				    toState+toModelDelta+toStrandDelta,
				    newP);
	}
      if(toModel.transExists(0,0))
	{
	  cerr << "MODEL "<<toModelId<<" HAS 0->0 TRANSITION"<<endl;
	  throw "transitivity not implemented for intergenic models!";
	}
    }
}



HMM *Application::loadHMM(const StringPlus &filename)
{
  // Load the forward model and its mirror image
  HMM &forwardModel=*new HMM(filename);
  HMM &reverseModel=*new HMM(forwardModel);
  reverseModel.reverseComplement();
  int numStates=forwardModel.getNumStates();
  strandDelta=numStates-1;
  numSymbols=forwardModel.getAlphabet().getNumElements();

  // Combine the forward and reverse models into the same model
  HMM &compositeModel=*new HMM(2*numStates-1,0,forwardModel.getAlphabet());
  copyIn(forwardModel,compositeModel,0);
  copyIn(reverseModel,compositeModel,strandDelta);

  // Add transitions from the zero state to both models (equiprobable)
  for(int i=0 ; i<numStates ; ++i)
    {
      double p=forwardModel.getTransProb(0,i);
      double halfP=p/2;
      compositeModel.setTransProb(0,i,halfP);
      compositeModel.setTransProb(0,i+strandDelta,halfP);
    }

  // Add these transitions:
  //   forward[downstream-intergenic] -> reverse[downstream-intergenic]
  //   reverse[upstream-intergenic] -> forward[upstream-intergenic]
  //cerr << "linking models, using strandDelta=" << strandDelta << endl;//###
  linkModels(downstreamModelId,0,downstreamModelId,strandDelta,
	     compositeModel,0.5);
  linkModels(upstreamModelId,strandDelta,upstreamModelId,0,
	     compositeModel,0.5);
  
  // Clean up 
  delete &forwardModel;
  delete &reverseModel;
  return &compositeModel;
}



ostream &operator<<(ostream &os,ExonTypes &exonType)
{
  switch(exonType)
    {
    case INITIAL_EXON:   os << "initial-exon";   break;
    case INTERNAL_EXON:  os << "internal-exon";  break;
    case FINAL_EXON:     os << "final-exon";     break;
    case SINGLE_EXON:    os << "single-exon";    break;
    }
  return os;
}




void Exon::printOn(ostream &os)
{
  os << axisId << "\tunveil\t" << exonType << "\t"
     << begin+1 << "\t" << end << "\t.\t" << strand << "\t"
     << frame <<"\t"<< transcriptId << endl;
}



ostream &operator<<(ostream &os,Exon &exon)
{
  exon.printOn(os);
  return os;
}



void Application::adjustExonTypes()
{
}



void Application::produceOutput()
{
  int n=exons.size();
  for(int i=0 ; i<n ; ++i)
    {
      Exon *exon=exons[i];
      cout << *exon;
    }
}



double Application::getPathScore(vector<int> &path,Sequence &seq,HMM &hmm)
{
  //cout << "pathlen=" << path.size() << ", seqlen=" << seq.getLength() << ", first state=" << path[0] << ", last state=" << path[path.size()-1] << endl;
  double logP=0;
  int pathLength=path.size();
  for(int i=0 ; i<pathLength ; ++i)
    {
      int state=path[i];
      Symbol symbol=seq[i];
      logP+=log(hmm.getEmitProb(state,symbol));
      //cout << "\t"<<state<<"="<<hmm.getEmitProb(state,symbol);
      if(i+1<pathLength)
	{
	  logP+=log(hmm.getTransProb(state,path[i+1]));
	  //cout << "\t" << state<<"->"<<path[i+1]<<" : " << hmm.getTransProb(state,path[i+1]) << endl;
	}
    }
  return logP;
}



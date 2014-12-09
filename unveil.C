/****************************************************************
 unveil.C
 An HMM genefinder based on Salzberg's VEIL model.

 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#include <string>
#include <iostream>
#include <math.h>
#include <fstream>
#include "tigr++/TigrCommandLine.H"
#include "tigr++/TigrSet.H"
#include "tigr++/TigrRegex.H"
#include "tigr++/TigrVector.H"
#include "tigr++/TigrMap.H"
#include "FastViterbi.H"
#include "tigr++/TigrFastaReader.H"
#include "tigr++/TigrFile.H"
#include "tigr++/DnaAlphabet.H"
using namespace std;


enum ExonType
  {
    INITIAL_EXON,
    INTERNAL_EXON,
    FINAL_EXON,
    SINGLE_EXON
  };
ostream &operator<<(ostream &,ExonType &);


class Exon
{
  int begin, end, frame, transcriptId;
  ExonType exonType;
  TigrString axisId;
  char strand;
public:
  Exon(int begin,int end,int frame,int transcriptId,
       ExonType exonType,TigrString axisId,char strand)
    : begin(begin), end(end), frame(frame), transcriptId(transcriptId),
      exonType(exonType), axisId(axisId), strand(strand) {}
  void printOn(ostream &);
  int getTranscriptId() const {return transcriptId;}
  ExonType getExonType() const {return exonType;}
  void setExonType(ExonType t) {exonType=t;}
  char getStrand() const {return strand;}
};
ostream &operator<<(ostream &,Exon &);


class Application
{
  DnaAlphabet alphabet;
  bool allowPartials;
  TigrSet<int> startCodonStates, stopCodonStates, exonStates,
    upstreamIntergenicStates, downstreamIntergenicStates, frameshiftStates;
  TigrRegex deflineRegex, donorRegex, acceptorRegex;
  vector<HiddenMarkovModel*> submodels;
  int upstreamModelId, downstreamModelId, strandDelta, numSymbols;
  vector<int> submodelOrigins; // state ids in composite HMM
  int numSubmodels;
  TigrVector<Exon*> exons;
  TigrMap<int,int> exonsPerGene;

  void copyIn(HiddenMarkovModel &submodel,HiddenMarkovModel &compositeHMM,
	      int delta);
  void addToSet(TigrSet<int> &,int fromState,int numStates);
  void loadSubmodels(const TigrString &filename);
  HiddenMarkovModel *loadHMM(const TigrString &);
  void linkOppStrandModels(int fromModelId,int fromDelta,char fromStrand,
			   int toModelId,int toDelta,char toStrand,
			   HiddenMarkovModel &compositeHMM,double metaP);
  void adjustExonTypes();
  void produceOutput();
  /*void closeTransitivity(int fromState,int toModelId,int toStrandDelta,
    double fromP,HiddenMarkovModel &compositeHMM);*/
  double getPathScore(vector<int> &path,Sequence &,HiddenMarkovModel &);
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
    TigrCommandLine cmd(argc,argv,"p");
    if(cmd.numArgs()!=3) 
      throw string("unveil [-p] <*.hmm> <submodels.txt> <*.fasta>\n"
"\t-p = allow partial genes");
    TigrString hmmFilename=cmd.arg(0);
    TigrString submodelsFilename=cmd.arg(1);
    TigrString fastaFilename=cmd.arg(2);
    allowPartials=cmd.option('p');

    time_t t=time(NULL);
    loadSubmodels(submodelsFilename);
    HiddenMarkovModel &hmm=*loadHMM(hmmFilename);
    FastViterbi viterbi(hmm);

    TigrFastaReader fastaReader(fastaFilename);
    TigrString defline, sequence;
    int transcriptId=0;
    while(fastaReader.nextSequence(defline,sequence))
      {
	if(!deflineRegex.search(defline))
	  throw TigrString("bad defline: ")+defline;
	TigrString axisId=deflineRegex[1];
	Sequence seq(sequence,alphabet);
	vector<int> &path=
	  *viterbi.getPath_Unveil(seq,frameshiftStates,exonStates,
				  startCodonStates,strandDelta);
	int len=path.size();
	int exonBegin=-1, frame=-9999999;
	bool inExon=false;
	ExonType exonType=INTERNAL_EXON;
	for(int i=0 ; i<len ; ++i)
	  {
	    int state=path[i];

	    // Detect forward-strand genes:

	    if(state<=strandDelta)
	      {
		if(!inExon && startCodonStates.isMember(state))
		  {
		    ++transcriptId;
		    exonsPerGene[transcriptId]=0;
		    cout << axisId << "\tunveil\tstart-codon\t"
			 << i+1 << "\t" << i+3 << "\t.\t+\t.\ttransgrp="
			 << transcriptId <<";"<< endl;
		    exonBegin=i;
		    inExon=true;
		    frame=0;
		    exonType=INITIAL_EXON;
		  }
		else if(!inExon && exonStates.isMember(state))
		  {
		    inExon=true;
		    exonBegin=i;
		    exonType=INTERNAL_EXON;

		    // ### This hack depends on knowledge of the loaded model!
		    if(i+1<len)
		      if(path[i+1]==116)
			{
			  if(i+2<len)
			    if(path[i+2]==117) frame=0;
			    else frame=1;
			}
		      else frame=2;
		    // ### end of hack
		  }
		else if(inExon && !exonStates.isMember(state))
		  {
		    if(i>0 && stopCodonStates.isMember(path[i-1]))
		      exonType=(exonType==INITIAL_EXON?SINGLE_EXON:FINAL_EXON);
		    exons.push_back(new Exon(exonBegin,i,frame,transcriptId,
					     exonType,axisId,'+'));
		    inExon=false;
		    exonType=INTERNAL_EXON;
		    ++exonsPerGene[transcriptId];
		  }
	      }

	    // Detect reverse-strand genes:

	    else
	      {
		if(!inExon && stopCodonStates.isMember(state-strandDelta))
		  {
		    ++transcriptId;
		    exonsPerGene[transcriptId]=0;
		    cout << axisId << "\tunveil\tstop-codon\t"
			 << i+1 << "\t" << i+3 << "\t.\t-\t.\ttransgrp="
			 << transcriptId <<";"<< endl;
		    exonBegin=i;
		    inExon=true;
		    exonType=FINAL_EXON;
		  }
		else if(!inExon && exonStates.isMember(state-strandDelta))
		  {
		    inExon=true;
		    exonBegin=i;
		    exonType=INTERNAL_EXON;
		  }
		else if(inExon && !exonStates.isMember(state-strandDelta))
		  {
		    // ### This hack depends on knowledge of the loaded model!
		    if(state-strandDelta==10) frame=0;
		    else
		      if(i-2>=0)
			if(path[i-2]-strandDelta==116)
			  {
			    if(i-3>=0)
			      if(path[i-3]-strandDelta==117) frame=0;
			      else frame=1;
			  }
			else frame=2;
		    // ### end of hack

		    if(i>0 && startCodonStates.isMember(path[i-1]-strandDelta))
		      exonType=(exonType==FINAL_EXON?SINGLE_EXON:INITIAL_EXON);
		    exons.push_back(new Exon(exonBegin,i,frame,transcriptId,
					     exonType,axisId,'-'));
		    inExon=false;
		    exonType=INTERNAL_EXON;
		    ++exonsPerGene[transcriptId];
		  }
	      }
	  }
	cout << "#score="<<getPathScore(path,seq,hmm)<<endl;
	delete &path;
      }

    adjustExonTypes();

    cout << "##gff-version 2\n"
	 << "##date " << ctime(&t)
	 << "##source-version unveil 1.0 (May 2003)\n"
	 << "##contact bmajoros@tigr.org\n"
	 << "##type DNA\n";
    produceOutput();
    cout << "# For information about this file format, see:" << endl
	 << "#   http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml" 
	 << endl;

    return 0;
  }



void Application::loadSubmodels(const TigrString &filename)
{
  submodelOrigins.push_back(0);
  int nextFreeState=1;
  TigrFile file(filename);
  numSubmodels=0;
  while(!file.eof())
    {
      TigrString line=file.readLine();
      if(file.eof()) break;
      vector<TigrString> &fields=*line.getFields("= \t\n[]");
      int numFields=fields.size();
      if(numFields==1)
	throw TigrString("Syntax error in submodel list file: ")+line;
      if(numFields>0)
	{
	  int modelId=fields[0].asInt();
	  TigrString secondField=fields[1];
	  TigrString modelFilename=fields[numFields-1];
	  HiddenMarkovModel *submodel=
	    new HiddenMarkovModel(modelFilename,alphabet);
	  while(modelId>=submodels.size()) submodels.push_back(NULL);
	  submodels[modelId]=submodel;
	  int submodelStates=submodel->countStates()-1;
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



void Application::addToSet(TigrSet<int> &theSet,int fromState,int numStates)
{
  int end=fromState+numStates;
  for(int i=fromState ; i<end ; ++i)
    theSet.insert(i);
}



void Application::copyIn(HiddenMarkovModel &submodel,
			 HiddenMarkovModel &compositeHMM,int delta)
{
  // Copy transitions
  int numStates=submodel.countStates();
  for(int i=1 ; i<numStates ; ++i)
    {
      for(int j=1 ; j<numStates ; ++j)
	compositeHMM.setTransitionProb(i+delta,j+delta,
				  submodel.getTransitionProb(i,j));
      compositeHMM.setTransitionProb(i+delta,0,
				submodel.getTransitionProb(i,0));
      compositeHMM.setTransitionProb(0,i+delta,
				submodel.getTransitionProb(0,i));
    }
  
  // Copy emissions
  for(int i=1 ; i<numStates ; ++i)
    for(int j=0 ; j<numSymbols ; ++j)
      compositeHMM.setEmissionProb(i+delta,Symbol(j),
			       submodel.getEmissionProb(i,Symbol(j)));
}



void Application::linkOppStrandModels(int fromModelId,int fromStrandDelta,
				      char fromStrand,int toModelId,
				      int toStrandDelta,char toStrand,
				      HiddenMarkovModel &compositeHMM,
				      double metaP)
{
  /*
    This method takes every state in the "from" model having a
    transition to zero and redirects that transition so that
    it points to the state(s) that the zero in the "to" model
    can transition to.
   */

  HiddenMarkovModel &fromModel=*submodels[fromModelId];
  HiddenMarkovModel &toModel=*submodels[toModelId];
  int fromModelDelta=submodelOrigins[fromModelId]-1;
  int toModelDelta=submodelOrigins[toModelId]-1;
  int numToModelStates=toModel.countStates();
  int numFromModelStates=fromModel.countStates();
  int numAllStates=compositeHMM.countStates();
  int skipBegin=fromModelDelta+fromStrandDelta+1;
  int skipEnd=fromModelDelta+fromStrandDelta+numFromModelStates;
  double intoZeroSum=0;
  if(toStrand=='-')
    for(int toState=1 ; toState<numToModelStates ; ++toState)
      intoZeroSum+=
	toModel.getTransitionProb(toState,0);

  for(int fromState=1 ; fromState<numFromModelStates ; ++fromState)
    {
      int adjustedFrom=fromState+fromModelDelta+fromStrandDelta;
      switch(fromStrand)
	{
	case '+':
	  if(!fromModel.doesTransitionExist(fromState,0)) continue;
	  break;
	case '-':
	  if(!fromModel.doesTransitionExist(0,fromState)) continue;
	  break;
	}

      double residue=0;
      for(int i=1 ; i<numAllStates ; ++i)
	{
	  if(i>=skipBegin && i<skipEnd) continue;
	  double halfOldP=compositeHMM.getTransitionProb(adjustedFrom,i)/2;
	  compositeHMM.setTransitionProb(adjustedFrom,i,halfOldP);
	  residue+=halfOldP;
	}
      for(int toState=1; toState<numToModelStates ; ++toState)
	{
	  double toP=(toStrand=='+' ?
		      toModel.getTransitionProb(0,toState) :
		      toModel.getTransitionProb(toState,0)/intoZeroSum);
	  if(toP==0) continue;
	  double newP=residue*toP;
	  if(compositeHMM.doesTransitionExist(fromState+fromModelDelta+
					 fromStrandDelta,
					 toState+toModelDelta+toStrandDelta))
	    throw TigrString("ERROR: transition already exists!");
	  compositeHMM.setTransitionProb(fromState+fromModelDelta+
					 fromStrandDelta,
					 toState+toModelDelta+toStrandDelta,
					 newP);
	}
      if(toModel.doesTransitionExist(0,0))
	{
	  cerr << "MODEL "<<toModelId<<" HAS 0->0 TRANSITION"<<endl;
	  throw "transitivity not implemented for intergenic models!";
	}
    }

}



/*
void Application::closeTransitivity(int fromState,int intermediateModelId,
				    int toStrandDelta,double fromP,
				    HiddenMarkovModel &compositeHMM)
{
  cout << "closing transitivity: "<<fromState<<"->"<<intermediateModelId
       << "("<<fromP<<")"<<endl;
  for(int toModelId=1 ; toModelId<=numSubmodels ; ++toModelId)
    if(metaModel->doesTransitionExist(intermediateModelId,toModelId))
      {
	double metaP=
	  metaModel->getTransitionProb(intermediateModelId,toModelId);
	HiddenMarkovModel &toModel=*submodels[toModelId];
	int toModelDelta=submodelOrigins[toModelId]-1;
	int numToModelStates=toModel.countStates();
	for(int toState=1 ; toState<numToModelStates ; ++toState)
	  {
	    if(!toModel.doesTransitionExist(0,toState)) continue;
	    double toP=toModel.getTransitionProb(0,toState);
	    double newP=fromP*metaP*toP;
	    if(compositeHMM.doesTransitionExist(fromState,
					toState+toModelDelta+toStrandDelta))
	      throw TigrString("ERROR: transition already exists!");
	    compositeHMM.setTransitionProb(fromState,
				      toState+toModelDelta+toStrandDelta,
				      newP);
	  }
	if(toModel.doesTransitionExist(0,0))
	  {
	    cout << "MODEL "<<toModelId<<" HAS 0->0 TRANSITION"<<endl;
	    double newP=fromP*metaP*toModel.getTransitionProb(0,0);
	    closeTransitivity(fromState,toModelId,toStrandDelta,
			      newP,compositeHMM,submodels);
	  }
      }
}
*/



HiddenMarkovModel *Application::loadHMM(const TigrString &filename)
{
  // Load the forward model and its mirror image
  HiddenMarkovModel &forwardModel=*new HiddenMarkovModel(filename,alphabet);
  HiddenMarkovModel &reverseModel=*new HiddenMarkovModel(forwardModel);
  reverseModel.reverseComp();

  int numStates=forwardModel.countStates();
  strandDelta=numStates-1;
  numSymbols=alphabet.getNumElements();

  // Combine the forward and reverse models into the same model
  HiddenMarkovModel &compositeModel=
    *new HiddenMarkovModel(alphabet,2*numStates-1);
  copyIn(forwardModel,compositeModel,0);
  copyIn(reverseModel,compositeModel,strandDelta);
  compositeModel.setTransitionProb(142,0,
     compositeModel.getTransitionProb(20,0)/2);//### HACK!!!
  compositeModel.setTransitionProb(161,160,
     compositeModel.getTransitionProb(1,2));//### HACK!!!
  compositeModel.setTransitionProb(161,161,
     compositeModel.getTransitionProb(1,1));//### HACK!!!

  /* Add these transitions:
     forward[downstream-intergenic] -> reverse[downstream-intergenic]
     reverse[upstream-intergenic] -> forward[upstream-intergenic]
  */
  linkOppStrandModels(downstreamModelId,0,'+',downstreamModelId,strandDelta,
		      '-',compositeModel,0.5);
  linkOppStrandModels(upstreamModelId,strandDelta,'-',upstreamModelId,0,'+',
		      compositeModel,0.5);

  // If partial genes are to be allowed, link zero state to all states
  numStates=compositeModel.countStates();
  if(allowPartials)
    {
      double p=1.0/(numStates-1);
      for(int i=1 ; i<numStates ; ++i)
	compositeModel.setTransitionProb(0,i,p);
      p=0.00001;
      double q=1-p;
      for(int i=2 ; i<numStates ; ++i)
	if(compositeModel.getTransitionProb(i,1)==0.0)
	{
	for(int j=0 ; j<numStates ; ++j)
	  compositeModel.setTransitionProb(i,j,
	    q*compositeModel.getTransitionProb(i,j));
	compositeModel.setTransitionProb(i,0,p);
	}
    }

  compositeModel.setTransitionProb(142,1,
    compositeModel.getTransitionProb(20,1)); // ### HACK!

  // Clean up 
  delete &forwardModel;
  delete &reverseModel;
  return &compositeModel;
}



ostream &operator<<(ostream &os,ExonType &exonType)
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
     << frame <<"\ttransgrp="<< transcriptId << ";" << endl;
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




double Application::getPathScore(vector<int> &path,Sequence &seq,
				 HiddenMarkovModel &hmm)
{
  double logP=0;
  int pathLength=path.size();
  for(int i=0 ; i<pathLength ; ++i)
    {
      int state=path[i];
      Symbol symbol=seq[i];
      logP+=log(hmm.getEmissionProb(state,symbol));
      if(i+1<pathLength) 
	logP+=log(hmm.getTransitionProb(state,path[i+1]));
    }
  return logP;
}


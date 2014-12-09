/**************************************************************
baum-welch.C : For training an HMM using Expectation Maximization
bmajoros@tigr.org 1/1/2003

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/

#include <string>
#include <iostream>
#include <fstream>
#include "tigr++/TigrFastaReader.H"
#include "tigr++/Sequence.H"
#include "tigr++/TigrRandom.H"
#include "tigr++/DnaAlphabet.H"
#include "tigr++/TigrCommandLine.H"
#include "tigr++/TigrVector.H"
#include "HMMreader.H"
#include "HMMbuilder.H"
#include "BaumWelch.H"
using namespace std;

class Application
{
  TigrVector<Sequence*> *trainingSet;
  DnaAlphabet alphabet;

  TigrVector<Sequence*> *load(TigrString);
public:
  int go(int argc,char *argv[]);
};


int main(int argc,char *argv[])
{
  try
    {
      Application app;
      return app.go(argc,argv);
    }
  catch(const char *msg)
    {
      cerr << msg << endl;
      return -1;
    }
  catch(const TigrString &s)
    {
      cerr << s.c_str() << endl;
      return -1;
    }
  catch(...)
    {
      cerr << "unknown exception caught in main()" << endl;
      return -1;
    }
  return 0;
}



void usage()
{
  throw "baum-welch <structure.hmms> <examples.fasta> <#iterations> <outfile>";
}


int Application::go(int argc,char *argv[])
{
  // Process command line
  TigrCommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=4) usage();
  TigrString structureFile=cmd.arg(0);
  TigrString trainFilename=cmd.arg(1);
  int maxIterations=cmd.arg(2).asInt();
  TigrString outfile=cmd.arg(3);
  randomize();

  // Load HMM structure
  HMMreader reader(alphabet);
  HiddenMarkovModel *hmm=reader.read(structureFile);
  //cerr << "Base model:\n" << *hmm << endl;

  // Load the training set
  trainingSet=load(trainFilename);

  // Train the HMM using the Baum-Welch algorithm
  cerr << "Training..." << endl;
  BaumWelch(*hmm,maxIterations,*trainingSet);
  cerr << "Final model:\n" << *hmm << endl;

  // Save the results
  hmm->save(outfile);

  return 0;
}


TigrVector<Sequence*> *Application::load(TigrString filename)
{
  TigrVector<Sequence*> *v=new TigrVector<Sequence*>;

  TigrFastaReader reader(filename);
  TigrString defline, sequence;
  while(reader.nextSequence(defline,sequence))
    {
      Sequence *seq=new Sequence(sequence,alphabet);
      v->push_back(seq);
    }
  return v;
}







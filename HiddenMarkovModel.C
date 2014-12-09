/**************************************************************
HiddenMarkovModel.C
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#include <fstream>
#include "HiddenMarkovModel.H"
#include "tigr++/TigrRandom.H"
#include <iostream>
#include "FastForward.H"
#include "HMMGraph.H"
using namespace std;


HiddenMarkovModel::HiddenMarkovModel(Alphabet &alphabet,int numStates)
  : numStates(numStates), emissionProb(numStates,alphabet.getNumElements()),
    transitionProb(numStates,numStates), alphabet(alphabet)
{
  emissionProb.setAllTo(0.0);
  transitionProb.setAllTo(0.0);
}



HiddenMarkovModel::HiddenMarkovModel(const HiddenMarkovModel &other)
  : emissionProb(other.emissionProb), numStates(other.numStates), 
    transitionProb(other.transitionProb), alphabet(other.alphabet)
{
}



HiddenMarkovModel::HiddenMarkovModel(const TigrString &filename,
				     Alphabet &alphabet)
  : emissionProb(0,0), transitionProb(0,0), alphabet(alphabet)
{
  load(filename);
}



HiddenMarkovModel::HiddenMarkovModel(istream &is,Alphabet &alphabet)
  : transitionProb(0,0), emissionProb(0,0), alphabet(alphabet)
{
  load(is);
}



void HiddenMarkovModel::printOn(ostream &os)
{
  os << "Transitions:\n";
  for(int i=0 ; i<numStates ; ++i)
    {
      for(int j=0 ; j<numStates ; ++j)
	os << '\t' << transitionProb[i][j];
      os << '\n';
    }
  int numSymbols=alphabet.getNumElements();
  os << "\nEmissions:\n";
  for(int i=0 ; i<numStates ; ++i)
    {
      for(int j=0 ; j<numSymbols ; ++j)
	os << '\t' << emissionProb[i][j];
      os << '\n';
    }  
  os << endl;
}



ostream &operator<<(ostream &os,HiddenMarkovModel &model)
{
  model.printOn(os);
  return os;
}



double HiddenMarkovModel::getTransitionProb(int from,int to)
{
  return transitionProb[from][to];
}



void HiddenMarkovModel::load(const TigrString &fname)
{
  ifstream is(fname.c_str());
  if(!is.good()) throw TigrString("Can't open ")+fname;
  load(is);
}



void HiddenMarkovModel::load(istream &is)
{
  int numSymbols=alphabet.getNumElements();

  // Load number of states and allocate arrays of that size
  is >> numStates;
  emissionProb.resize(numStates,numSymbols);
  transitionProb.resize(numStates,numStates);

  // Read transition  && emission probabilities
  is >> transitionProb >> emissionProb;
}



void HiddenMarkovModel::normalizeTransitions()
{
  double sum;
  for(int i=0 ; i<numStates ; ++i)
    {
      sum=0;
      for(int j=0 ; j<numStates ; ++j)
	sum+=transitionProb[i][j];
      if(sum>0)
	for(int j=0 ; j<numStates ; ++j)
	  transitionProb[i][j]/=sum;
    }
}



bool HiddenMarkovModel::save(const TigrString &fname)
{
  ofstream os(fname.c_str());
  if(!os.good()) throw TigrString("Can't create ")+fname;
  bool success=save(os);
  return success;
}



bool HiddenMarkovModel::save(ostream &os)
{
  os << numStates << '\n' 
     << transitionProb << '\n' 
     << emissionProb << endl;
  return true;
}



void HiddenMarkovModel::reverseComp()
{
  double t;
  Symbol C=alphabet.lookup('C');
  Symbol T=alphabet.lookup('T');
  Symbol A=alphabet.lookup('A');
  Symbol G=alphabet.lookup('G');
  swap(G,C);
  swap(A,T);
  normalizeEmissions();
  for(int i=0 ; i<numStates ; ++i)
    for(int j=i+1 ; j<numStates ; ++j)
      {
	t=transitionProb[j][i];
	transitionProb[j][i]=transitionProb[i][j];
	transitionProb[i][j]=t;
      }

  // NOTE:  do NOT call normalizeTransitions() here unless
  //        you really want to ... otherwise, the reverseComp
  //        operation will not be "reversible" (i.e., a
  //        symmetric binary relation).
}



bool HiddenMarkovModel::doesTransitionExist(int from,int to)
{
  bool exists=transitionProb[from][to]>0;
  return exists;
}



void HiddenMarkovModel::normalizeEmissions()
{
  double sum;
  int numSymbols=alphabet.getNumElements();
  for(int from=0 ; from<numStates ; ++from)
    {
      sum=0;
      for(int s=0 ; s<numSymbols ; ++s)
	sum+=emissionProb[from][s];
      if(sum>0)
	for(int s=0 ; s<numSymbols ; ++s)
	  emissionProb[from][s]/=sum;
    }
}



void HiddenMarkovModel::addTransition(int from,int to)
{
  setTransitionProb(from,to,Random0to1());
}



int HiddenMarkovModel::countStates()
{
  return numStates;
}



void HiddenMarkovModel::swap(Symbol s1,Symbol s2)
{
  double t;
  for(int state=0 ; state<numStates ; ++state)
    {
      t=emissionProb[state][s1];
      emissionProb[state][s1]=emissionProb[state][s2];
      emissionProb[state][s2]=t;
    }
}



void HiddenMarkovModel::setTransitionProb(int from,int to,double prob)
{
  transitionProb[from][to]=prob;
}



double HiddenMarkovModel::getEmissionProb(int state,Symbol s)
{
  return emissionProb[state][s];
}



Alphabet &HiddenMarkovModel::getAlphabet()
{
  return alphabet;
}



double HiddenMarkovModel::getLogP(Sequence &seq,int begin,int len)
{
  if(len<0) len=seq.getLength();
  HMMGraph hmmGraph(*this);
  FastForward f(*this,hmmGraph,seq,begin,len);
  return f.getLogP();
}



void HiddenMarkovModel::initEmissionProbs()
{
  double sum;
  int numSymbols=alphabet.getNumElements();
  for(int i=0 ; i<numStates ; ++i)
    {
      sum=0.0;
      for(int j=0 ; j<numSymbols ; ++j)
	{
	  double r=RandomFloat(0.20,0.80);
	  sum+=r;
	  emissionProb[i][j]=r;
	}
      if(sum>0)
	for(int j=0 ; j<numSymbols ; ++j)
	  { emissionProb[i][j]/=sum; }
    }
}



void HiddenMarkovModel::setEmissionProb(int state,Symbol s,double prob)
{
  emissionProb[state][s]=prob;
}




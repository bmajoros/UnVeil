/**************************************************************
FastViterbi.C
bmajoros@tigr.org

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#include <math.h>
#include <float.h>
#include <iostream>
#include "FastViterbi.H"
#include "tigr++/Constants.H"
using namespace std;

#define MOST_NEGATIVE -DBL_MAX

inline double safeAdd(double a,double b)
{
  if(a==MOST_NEGATIVE || b==MOST_NEGATIVE) return MOST_NEGATIVE;
  else return a+b;
}


FastViterbi::FastViterbi(HiddenMarkovModel &hmm)
  : numStates(hmm.countStates()), hmmGraph(hmm)
{
  // ctor
}



TigrVector<int> *FastViterbi::buildPath(int L,int finalState,
					TigrArray2D<short> &linkBack)
{
  int currentState=finalState;
  TigrVector<int> &path=*new TigrVector<int>(L);
  path[L-1]=finalState;
  for(int pos=L-1 ; pos>0 ; --pos)
    {
      currentState=linkBack[currentState][pos+1];
      path[pos-1]=currentState;
    }
  return &path;
}



/*
  This method returns the state sequence of
  the most probable path (not including the
  start state).
 */
TigrVector<int> *FastViterbi::getPath(Sequence &sequence)
{
  HiddenMarkovModel &hmm=hmmGraph.getHMM();
  int L=sequence.getLength(), Lplus1=L+1;
  TigrArray2D<short> linkBack(numStates,Lplus1);
  TigrArray2D<double> m(numStates,Lplus1);
  m.setAllTo(NEGATIVE_INFINITY);
  m[0][0]=0.0;
  for(int i=1 ; i<=L ; ++i)
    {
      Symbol s=sequence[i-1];
      TigrVector<StateDoublePair> &statesEmittingS=
	hmmGraph.statesEmitting(s);
      int numEmitting=statesEmittingS.size();
      for(int l=0 ; l<numEmitting ; ++l)
	{
	  StateDoublePair &currentPair=statesEmittingS[l];
	  TigrVector<StateDoublePair> &precedingStates=
	    hmmGraph.statesPreceding(currentPair.state);
	  int numPreceding=precedingStates.size();
	  double bestP=log(0.0);
	  int bestPredecessor;
	  for(int j=0 ; j<numPreceding ; ++j)
	    {
	      StateDoublePair precedingPair=precedingStates[j];
	      double inductiveP=m[precedingPair.state][i-1];
	      double newP=safeAdd(inductiveP,precedingPair.logP);
	      if(newP>bestP)
		{
		  bestP=newP;
		  bestPredecessor=precedingPair.state;
		}
	    }
	  m[currentPair.state][i]=safeAdd(bestP,currentPair.logP);
	  linkBack[currentPair.state][i]=bestPredecessor;
	}
    }
  return buildPath(L,0,linkBack);
}



TigrVector<int> *FastViterbi::getPath_Unveil(Sequence &sequence,
					     TigrSet<int> &frameshiftStates,
					     TigrSet<int> &exonStates,
					     TigrSet<int> &startCodonStates,
					     int strandDelta)
{
  HiddenMarkovModel &hmm=hmmGraph.getHMM();
  const int L=sequence.getLength();
  TigrArray2D<short> linkBack(numStates,L+1);
  TigrArray2D<double> m(numStates,L+1);
  TigrArray2D<short> frameshifts(numStates,L+1);
  m.setAllTo(NEGATIVE_INFINITY);
  m[0][0]=0.0;
  frameshifts.setAllTo(0); // ### could be made faster
  for(int i=1 ; i<=L ; ++i)
    {
      Symbol s=sequence[i-1];
      TigrVector<StateDoublePair> &statesEmittingS=
	hmmGraph.statesEmitting(s);
      int numEmitting=statesEmittingS.size();
      for(int l=0 ; l<numEmitting ; ++l)
	{
	  StateDoublePair &currentPair=statesEmittingS[l];
	  TigrVector<StateDoublePair> &precedingStates=
	    hmmGraph.statesPreceding(currentPair.state);
	  int numPreceding=precedingStates.size();
	  if(!numPreceding) continue;

	  // #### UNVEIL-specific code:
	  int adjustedCurrentState=currentPair.state;
	  if(adjustedCurrentState>=strandDelta)
	    adjustedCurrentState-=strandDelta;
           // ###

	  double bestP=log(0.0);
	  int bestPredecessor=-1;
	  for(int j=0 ; j<numPreceding ; ++j)
	    {
	      StateDoublePair precedingPair=precedingStates[j];

	      // #### UNVEIL-specific code:
	      int adjustedPrecedingState=precedingPair.state;
	      if(adjustedPrecedingState>=strandDelta)
		adjustedPrecedingState-=strandDelta;
	      // ### THIS HACK IS HIGHLY MODEL-DEPENDENT!
	      if(adjustedCurrentState>=30 && adjustedCurrentState<=113 &&
		 (adjustedPrecedingState<30 || adjustedPrecedingState>113) 
		 && frameshifts[precedingPair.state][i-1]%3)
		continue;
	      // ####

	      double inductiveP=m[precedingPair.state][i-1];
	      double newP=safeAdd(inductiveP,precedingPair.logP);
	      if(newP>bestP)
		{
		  bestP=newP;
		  bestPredecessor=precedingPair.state;
		}

	      
	    }
	  if(bestPredecessor<0) continue;//### 5/8/03

	  // #### UNVEIL-specific code:
	  frameshifts[currentPair.state][i]=
	    frameshifts[bestPredecessor][i-1];
	  if(frameshiftStates.isMember(adjustedCurrentState))
	    ++frameshifts[currentPair.state][i];
	  // ###

	  m[currentPair.state][i]=safeAdd(bestP,currentPair.logP);
	  linkBack[currentPair.state][i]=bestPredecessor;
	}
    }

  int bestState=1;
  double bestScore=NEGATIVE_INFINITY;
  for(int i=0 ; i<numStates ; ++i)
    if(m[i][L-1]>bestScore)
      {
	bestState=i;
	bestScore=m[i][L-1];
      }

  return buildPath(L,1,linkBack);
}




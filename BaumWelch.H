/****************************************************************
 BaumWelch.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/

#ifndef INCL_BaumWelch_H
#define INCL_BaumWelch_H

#include "ForwardAlgorithm.H"
#include "BackwardAlgorithm.H"
#include <math.h>
#include <iostream>
#include <strstream>
#include <fstream>
#include <float.h>
#include "tigr++/TigrVector.H"
#include "tigr++/TigrArray1D.H"
#include "tigr++/Sequence.H"
#include "tigr++/TigrProgress.H"
using namespace std;


class BaumWelch
{
public:
  BaumWelch(HiddenMarkovModel &,long maxIterations,
	    TigrVector<Sequence*> &trainingSet);
private:
  Alphabet &alphabet;
  TigrProgress progress;
  TigrArray1D<double> likelihoods, scalingRatios;
  TigrArray2D<double> A; // expected transition counts
  TigrArray2D<double> E; // expected emission counts
  HiddenMarkovModel &hmm;
  long maxIterations;
  double logLikelihoodThreshold, sum;
  TigrVector<Sequence*> &trainingSet;
  int numHmmStates, numTrain, alphabetSize;
  
  void reviseExpectedEmitCounts(int sequenceLength,
				ForwardAlgorithm &,BackwardAlgorithm &,
				Sequence &);
  void reviseExpectedTransCounts(int sequenceLength,ForwardAlgorithm &,
				 BackwardAlgorithm &,
				 Sequence &);
  void mainAlgorithm();
};


#endif

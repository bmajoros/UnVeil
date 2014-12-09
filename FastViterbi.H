/**************************************************************
FastViterbi.H
bmajoros@tigr.org

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_FastViterbi_H
#define INCL_FastViterbi_H

#include <iostream>
#include "tigr++/TigrSet.H"
#include "HiddenMarkovModel.H"
#include "tigr++/Sequence.H"
#include "tigr++/TigrBitSet.H"
#include "HMMGraph.H"
#include "tigr++/TigrVector.H"
#include "tigr++/TigrArray2D.H"
using namespace std;


class FastViterbi
{
public:
  FastViterbi(HiddenMarkovModel &);
  TigrVector<int> *getPath(Sequence &);
  TigrVector<int> *getPath_Unveil(Sequence &,
				  TigrSet<int> &frameshiftStates,
				  TigrSet<int> &exonStates,
				  TigrSet<int> &startCodonStates,
				  int strandDelta);
private:
  int numStates;
  HMMGraph hmmGraph;
  TigrVector<int> *buildPath(int L,int pathEndState,
			     TigrArray2D<short> &ptr);
};


#endif

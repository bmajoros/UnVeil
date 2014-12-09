/**************************************************************
BaumWelch.C : the Baum-Welch Expectation Maximization algorithm
bmajoros@tigr.org

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/

#include "BaumWelch.H"

BaumWelch::BaumWelch(HiddenMarkovModel &hmm,long maxIterations,
		     TigrVector<Sequence*> &trainingSet)
  : hmm(hmm), 
    numHmmStates(hmm.countStates()),alphabet(hmm.getAlphabet()),
    numTrain(trainingSet.size()),likelihoods(trainingSet.size()),
    trainingSet(trainingSet),A(hmm.countStates(),hmm.countStates()),
    alphabetSize(alphabet.getNumElements()),scalingRatios(1),
    maxIterations(maxIterations),
    E(hmm.countStates(),alphabet.getNumElements())
{
  mainAlgorithm();
}



/***************************************************************
  This is the main Baum-Welch Expectation Maximization algorithm
***************************************************************/
void BaumWelch::mainAlgorithm()
{
  double logLikelihood, deltaLikelihood=logLikelihoodThreshold+1, 
    oldLogLikelihood=0;

  /*************************************************************
    Various initialization steps
  *************************************************************/
  int numEmptyStrings=0;
  for(int j=0 ; j<numTrain ; ++j)
    {
      Sequence &sequence=*trainingSet[j];
      if(sequence.getLength()==0) ++numEmptyStrings;
    }
  bool allowEmptyStrings=hmm.doesTransitionExist(0,0);
  hmm.setTransitionProb(0,0,0);
  progress.start(maxIterations);

  /***********************************************************
    Main loop
  ***********************************************************/
  for(int i=0 ; i<maxIterations ; ++i)
    {
      A.setAllTo(0); // don't use pseudocounts--they cause likelihood drops!
      E.setAllTo(0); // (unless you re-do the math to account for them...)
      for(int symbol=0 ; symbol<alphabetSize ; ++symbol) E[0][symbol]=0;
      
      double likelihoodCorrection=0.0;
      for(int j=0 ; j<numTrain ; ++j)
	{
	  Sequence &sequence=*trainingSet[j];
	  int sequenceLength=sequence.getLength();
	  if(sequenceLength>0)
	    {
	      /***********************************************************
		Run the Forward algorithm to compute f(k,i)=P[a model M in
		the initial state will emit sequence x(1)..x(i) and reside
		in state k when emitting xi at time i]:
	       ***********************************************************/
	      ForwardAlgorithm f(hmm,sequence);
	      const TigrArray1D<double> &scalingValues=
		f.getScalingFactors();
	      likelihoods[j]=f.getScaledP();
	      for(int i=1 ; i<=sequenceLength ; ++i)
		likelihoodCorrection+=log(scalingValues[i]);

	      /***********************************************************
		Run the Backward algorithm to compute b(k,i)=P[M will 
		next emit the string x(i+1)..x(L) and then terminate |
		M is currently in state k]:
	       ***********************************************************/
	      BackwardAlgorithm b(f.getScalingFactors(),hmm,sequence);

	      /***********************************************************
		Compute scaling ratios so the scaling values can be
		properly cancelled during the updating of the A and E
		counts:
	       ***********************************************************/
	      scalingRatios.resize(sequenceLength+1);
	      scalingRatios[sequenceLength]=
		1/f.getScalingFactors()[sequenceLength];
	      for(int i=sequenceLength-1 ; i>0 ; --i)
		{
		  double r=scalingRatios[i+1] *
		    b.getScalingFactors()[i+1]/f.getScalingFactors()[i];
		  scalingRatios[i]=r;
		  if(isinf(r) || isnan(r) || r==0) 
		    cerr << "WARNING: scaling ratio "<<i<<"="
			 <<scalingRatios[i]; 
		}

	      /***********************************************************
		Update the expected number of transitions and emissions
		used in generating this sequence:
	       ***********************************************************/
	      reviseExpectedTransCounts(sequenceLength,f,b,sequence);
	      reviseExpectedEmitCounts(sequenceLength,f,b,sequence);
	    }
	}

      /********************************************************************
	Recompute the a's and e's from the A's and E's (i.e., probabilities 
	from counts)
      ********************************************************************/
      for(int k=1 ; k<numHmmStates ; ++k)
	{
	  sum=0;
	  for(Symbol s=0 ; s<alphabetSize ; ++s) sum+=E[k][s];
	  if(sum==0) throw TigrString("WARNING: HMM state ")+k+
		       " cannot emit any symbol; possibly unreachable.";
	  else for(Symbol s=0 ; s<alphabetSize ; ++s)
	    {
	      double newValue=E[k][s]/sum;
	      hmm.setEmissionProb(k,s,newValue);
	      if(newValue==1)
		cerr << "WARNING: hmm.getEmitProb(k,s)=1 for k=" 
		     << k << ", s=" << s << endl; 
	    }
	}
      for(int k=0 ; k<numHmmStates ; ++k)
	{
	  sum=0;
	  for(int ll=0 ; ll<numHmmStates ; ++ll) sum+=A[k][ll];
	  if(sum==0) for(int ll=0 ; ll<numHmmStates ; ++ll) 
	    hmm.setTransitionProb(k,ll,0);
	  else for(int l=0 ; l<numHmmStates ; ++l)
	    {
	      double p=A[k][l]/sum;
	      hmm.setTransitionProb(k,l,p);
	    }
	}

      /*********************************************************************
	Compute the log likelihood.  This is only necessary during 
        development and debugging, to ensure that recent changes to the 
        code have not changed the invariant that likelihood should increase
        monotonically (possibly except for tiny fluctuations due to 
        rounding in a digital computer):
      *********************************************************************/
      double logLikelihood=likelihoodCorrection;
      for(int j=0 ; j<numTrain ; ++j)
	if(trainingSet[j]->getLength()>0)
	  logLikelihood+=log(likelihoods[j]);
      deltaLikelihood=logLikelihood-oldLogLikelihood;
      oldLogLikelihood=logLikelihood;
      cerr << "log(likelihood)=" << logLikelihood;
      if(i>0)
	{
	  cerr << " " << progress.getProgress(i) << endl; 
	  if(deltaLikelihood<0)
	    cerr << "^--WARNING!  LIKELIHOOD DECREASED--^" << endl;
	}
      else cerr<<endl;
    }

  /*************************************************************************
   Allow start-state-self-transitions only if empty strings are to be 
   permitted
   ************************************************************************/
  if(allowEmptyStrings)
    {
      double p=double(numEmptyStrings)/double(numTrain);
      hmm.setTransitionProb(0,0,p);
      double q=1-p;
      for(int state=1 ; state<numHmmStates ; ++state)
	{
	  double rescaled=hmm.getTransitionProb(0,state)*q;
	  hmm.setTransitionProb(0,state,rescaled);
	}
    }
}



/*
  Update the expected emission counts for the HMM for a training sequence.
 */
void BaumWelch::reviseExpectedEmitCounts(int sequenceLength,
					 ForwardAlgorithm &f,
					 BackwardAlgorithm &b,
					 Sequence &sequence)
{
  double seqProb=f.getScaledP();
  for(int k=1 ; k<numHmmStates ; ++k)
    for(Symbol s=0 ; s<alphabetSize ; ++s)
      {
	sum=0.0;
	for(int i=1 ; i<=sequenceLength ; ++i)
	  if(s==sequence[i-1])
	    {
	      sum+=f.getScalingFactors()[i] * scalingRatios[i]
		* f(k,i) * b(k,i) / seqProb;
	      if(isinf(sum) || isnan(sum)) throw TigrString("sum=")+sum+
		  " in BaumWelch::reviseExpectedEmitCounts()";
	    }
	E[k][s]+=sum;
      }
}



/*
  Update the expected transition counts for the HMM for a training sequence.
 */
void BaumWelch::reviseExpectedTransCounts(int sequenceLength,
					  ForwardAlgorithm &f,
					  BackwardAlgorithm &b,
					  Sequence &sequence)
{
  double seqProb=f.getScaledP();
  for(int k=0 ; k<numHmmStates ; ++k)
    for(int l=1 ; l<numHmmStates ; ++l)
      {
	sum=0.0;
	for(int i=0 ; i<sequenceLength ; ++i)
	  {
	    sum+=f(k,i) * hmm.getTransitionProb(k,l) * 
	      hmm.getEmissionProb(l,sequence[i]) * b(l,i+1) / seqProb * 
	      scalingRatios[i+1];
	    if(isinf(sum) || isnan(sum)) throw TigrString("sum=")+sum+
		   " in BaumWelch::reviseExpectedTransCounts()";
	  }
	A[k][l]+=sum;
      }
  for(int k=1 ; k<numHmmStates ; ++k)
    A[k][0]+=f(k,sequenceLength) * hmm.getTransitionProb(k,0) / seqProb *
      scalingRatios[sequenceLength] * f.getScalingFactors()[sequenceLength];
  A[0][0]=0;
}




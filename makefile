CC 		= g++
LDFLAGS	        = -O
CFLAGS		= -O -w -fpermissive
OBJ		= obj
TIGRPLUSPLUS	= tigr++
UNVEIL		= .
LIBS		= 

all: \
	obj \
	deterministic-trainer \
	unveil \
	baum-welch \
	model-combiner \
	perl-files \
	warning \

.PHONY : perl-files
perl-files:
	chmod +x *.pl

.PHONY : warning
warning:
	@echo WARNING!  Do not attempt to train the gene-finder in the installation directory.  Read the included README file.  An auto-trainer script is provided for your convenience.

obj:
	mkdir obj

.PHONY : clean
clean:
	@rm obj/*.o


#---------------------------------------------------------
$(OBJ)/regex.o:\
		$(TIGRPLUSPLUS)/regex.h \
		$(TIGRPLUSPLUS)/regex.c
	gcc $(CFLAGS) -o $(OBJ)/regex.o -c \
		$(TIGRPLUSPLUS)/regex.c

$(OBJ)/TigrFastaReader.o: \
		$(TIGRPLUSPLUS)/TigrFastaReader.H \
		$(TIGRPLUSPLUS)/TigrFastaReader.C
	g++ $(CFLAGS) -o $(OBJ)/TigrFastaReader.o -c \
		$(TIGRPLUSPLUS)/TigrFastaReader.C

$(OBJ)/DnaAlphabet.o: \
		$(TIGRPLUSPLUS)/DnaAlphabet.H \
		$(TIGRPLUSPLUS)/DnaAlphabet.C
	g++ $(CFLAGS) -o $(OBJ)/DnaAlphabet.o -c \
		$(TIGRPLUSPLUS)/DnaAlphabet.C

$(OBJ)/model-combiner.o: \
		$(UNVEIL)/model-combiner.C \
		$(TIGRPLUSPLUS)/Alphabet.H \
		$(UNVEIL)/ForwardAlgorithm.H \
		$(UNVEIL)/BaumWelch.H \
		$(UNVEIL)/HiddenMarkovModel.H \
		$(TIGRPLUSPLUS)/Sequence.H
	g++ $(CFLAGS) -o $(OBJ)/model-combiner.o -c \
		$(UNVEIL)/model-combiner.C

$(OBJ)/deterministic-trainer.o: \
		$(UNVEIL)/deterministic-trainer.C \
		$(TIGRPLUSPLUS)/Alphabet.H \
		$(UNVEIL)/ForwardAlgorithm.H \
		$(UNVEIL)/BaumWelch.H \
		$(UNVEIL)/HiddenMarkovModel.H \
		$(TIGRPLUSPLUS)/Sequence.H
	g++ $(CFLAGS) -o $(OBJ)/deterministic-trainer.o -c \
		$(UNVEIL)/deterministic-trainer.C

$(OBJ)/baum-welch.o: \
		$(UNVEIL)/baum-welch.C \
		$(TIGRPLUSPLUS)/Alphabet.H \
		$(UNVEIL)/ForwardAlgorithm.H \
		$(UNVEIL)/BaumWelch.H \
		$(UNVEIL)/HiddenMarkovModel.H \
		$(TIGRPLUSPLUS)/Sequence.H
	g++ $(CFLAGS) -o $(OBJ)/baum-welch.o -c \
		$(UNVEIL)/baum-welch.C

$(OBJ)/Alphabet.o: \
		$(TIGRPLUSPLUS)/Alphabet.H \
		$(TIGRPLUSPLUS)/Alphabet.C
	g++ $(CFLAGS) -o $(OBJ)/Alphabet.o -c \
		$(TIGRPLUSPLUS)/Alphabet.C

$(OBJ)/HMMbuilder.o: \
		$(UNVEIL)/HMMbuilder.H \
		$(UNVEIL)/HMMbuilder.C
	g++ $(CFLAGS) -o $(OBJ)/HMMbuilder.o -c \
		$(UNVEIL)/HMMbuilder.C

$(OBJ)/HMMreader.o: \
		$(UNVEIL)/HMMreader.H \
		$(UNVEIL)/HMMreader.C
	g++ $(CFLAGS) -o $(OBJ)/HMMreader.o -c \
		$(UNVEIL)/HMMreader.C

$(OBJ)/BackwardAlgorithm.o: \
		$(UNVEIL)/HiddenMarkovModel.H \
		$(TIGRPLUSPLUS)/Sequence.H \
		$(UNVEIL)/BackwardAlgorithm.H \
		$(UNVEIL)/BackwardAlgorithm.C
	g++ $(CFLAGS) -o $(OBJ)/BackwardAlgorithm.o -c \
		$(UNVEIL)/BackwardAlgorithm.C

$(OBJ)/BaumWelch.o: \
		$(UNVEIL)/HiddenMarkovModel.H \
		$(TIGRPLUSPLUS)/Sequence.H \
		$(UNVEIL)/ForwardAlgorithm.H \
		$(UNVEIL)/BackwardAlgorithm.H \
		$(UNVEIL)/BaumWelch.H \
		$(UNVEIL)/BaumWelch.C
	g++ $(CFLAGS) -o $(OBJ)/BaumWelch.o -c \
		$(UNVEIL)/BaumWelch.C

$(OBJ)/ForwardAlgorithm.o: \
		$(TIGRPLUSPLUS)/Sequence.H \
		$(UNVEIL)/HiddenMarkovModel.H \
		$(UNVEIL)/ForwardAlgorithm.H \
		$(UNVEIL)/ForwardAlgorithm.C
	g++ $(CFLAGS) -o $(OBJ)/ForwardAlgorithm.o -c \
		$(UNVEIL)/ForwardAlgorithm.C

$(OBJ)/FastViterbi.o: \
		$(TIGRPLUSPLUS)/Sequence.H \
		$(UNVEIL)/HiddenMarkovModel.H \
		$(UNVEIL)/HMMGraph.H \
		$(UNVEIL)/FastViterbi.H \
		$(UNVEIL)/FastViterbi.C
	g++ $(CFLAGS) -o $(OBJ)/FastViterbi.o -c \
		$(UNVEIL)/FastViterbi.C

$(OBJ)/HiddenMarkovModel.o: \
		$(UNVEIL)/HiddenMarkovModel.H \
		$(UNVEIL)/HiddenMarkovModel.C \
		$(TIGRPLUSPLUS)/Alphabet.H \
		$(TIGRPLUSPLUS)/Symbol.H \
		$(TIGRPLUSPLUS)/TigrRandom.H
	g++ $(CFLAGS) -o $(OBJ)/HiddenMarkovModel.o -c \
		$(UNVEIL)/HiddenMarkovModel.C

$(OBJ)/Sequence.o: \
		$(TIGRPLUSPLUS)/Sequence.H \
		$(TIGRPLUSPLUS)/Symbol.H \
		$(TIGRPLUSPLUS)/Alphabet.H \
		$(TIGRPLUSPLUS)/Sequence.C
	g++ $(CFLAGS) -o $(OBJ)/Sequence.o -c \
		$(TIGRPLUSPLUS)/Sequence.C

$(OBJ)/Symbol.o: \
		$(TIGRPLUSPLUS)/Symbol.H \
		$(TIGRPLUSPLUS)/Symbol.C
	g++ $(CFLAGS) -o $(OBJ)/Symbol.o -c \
		$(TIGRPLUSPLUS)/Symbol.C

$(OBJ)/TigrRandom.o: \
		$(TIGRPLUSPLUS)/TigrRandom.H \
		$(TIGRPLUSPLUS)/TigrRandom.C
	g++ $(CFLAGS) -o $(OBJ)/TigrRandom.o -c \
		$(TIGRPLUSPLUS)/TigrRandom.C

$(OBJ)/TigrProgress.o: \
		$(TIGRPLUSPLUS)/TigrProgress.H \
		$(TIGRPLUSPLUS)/TigrProgress.C
	g++ $(CFLAGS) -o $(OBJ)/TigrProgress.o -c \
		$(TIGRPLUSPLUS)/TigrProgress.C

$(OBJ)/TigrCommandLine.o: \
		$(TIGRPLUSPLUS)/TigrCommandLine.H \
		$(TIGRPLUSPLUS)/TigrCommandLine.C
	g++ $(CFLAGS) -o $(OBJ)/TigrCommandLine.o -c \
		$(TIGRPLUSPLUS)/TigrCommandLine.C

$(OBJ)/TigrConfigFile.o: \
		$(TIGRPLUSPLUS)/TigrConfigFile.H \
		$(TIGRPLUSPLUS)/TigrConfigFile.C
	g++ $(CFLAGS) -o $(OBJ)/TigrConfigFile.o -c \
		$(TIGRPLUSPLUS)/TigrConfigFile.C

$(OBJ)/TigrStrTokenizer.o: \
		$(TIGRPLUSPLUS)/TigrStrTokenizer.H \
		$(TIGRPLUSPLUS)/TigrStrTokenizer.C
	g++ $(CFLAGS) -o $(OBJ)/TigrStrTokenizer.o -c \
		$(TIGRPLUSPLUS)/TigrStrTokenizer.C

$(OBJ)/TigrFile.o: \
		$(TIGRPLUSPLUS)/TigrFile.H \
		$(TIGRPLUSPLUS)/TigrFile.C
	g++ $(CFLAGS) -o $(OBJ)/TigrFile.o -c \
		$(TIGRPLUSPLUS)/TigrFile.C

$(OBJ)/TigrString.o: \
		$(TIGRPLUSPLUS)/TigrString.H \
		$(TIGRPLUSPLUS)/TigrString.C
	g++ $(CFLAGS) -o $(OBJ)/TigrString.o -c \
		$(TIGRPLUSPLUS)/TigrString.C

$(OBJ)/TigrRegex.o: \
		$(TIGRPLUSPLUS)/TigrRegex.H \
		$(TIGRPLUSPLUS)/TigrRegex.C
	g++ $(CFLAGS) -o $(OBJ)/TigrRegex.o -c \
		$(TIGRPLUSPLUS)/TigrRegex.C

$(OBJ)/TigrExceptions.o: \
		$(TIGRPLUSPLUS)/TigrExceptions.H \
		$(TIGRPLUSPLUS)/TigrExceptions.C
	g++ $(CFLAGS) -o $(OBJ)/TigrExceptions.o -c \
		$(TIGRPLUSPLUS)/TigrExceptions.C

$(OBJ)/TigrStacktrace.o: \
		$(TIGRPLUSPLUS)/TigrStacktrace.H \
		$(TIGRPLUSPLUS)/TigrStacktrace.C
	g++ $(CFLAGS) -o $(OBJ)/TigrStacktrace.o -c \
		$(TIGRPLUSPLUS)/TigrStacktrace.C

$(OBJ)/TigrBitSet.o: \
		$(TIGRPLUSPLUS)/TigrBitSet.H \
		$(TIGRPLUSPLUS)/TigrBitSet.C
	g++ $(CFLAGS) -o $(OBJ)/TigrBitSet.o -c \
		$(TIGRPLUSPLUS)/TigrBitSet.C

$(OBJ)/TigrGffReader.o: \
		$(TIGRPLUSPLUS)/TigrGffReader.H \
		$(TIGRPLUSPLUS)/TigrGffReader.C
	g++ $(CFLAGS) -o $(OBJ)/TigrGffReader.o -c \
		$(TIGRPLUSPLUS)/TigrGffReader.C
$(OBJ)/TigrGffExon.o: \
		$(TIGRPLUSPLUS)/TigrGffExon.H \
		$(TIGRPLUSPLUS)/TigrGffExon.C
	g++ $(CFLAGS) -o $(OBJ)/TigrGffExon.o -c \
		$(TIGRPLUSPLUS)/TigrGffExon.C
$(OBJ)/TigrGffTranscript.o: \
		$(TIGRPLUSPLUS)/TigrGffTranscript.H \
		$(TIGRPLUSPLUS)/TigrGffTranscript.C
	g++ $(CFLAGS) -o $(OBJ)/TigrGffTranscript.o -c \
		$(TIGRPLUSPLUS)/TigrGffTranscript.C
$(OBJ)/TigrGffFeature.o: \
		$(TIGRPLUSPLUS)/TigrGffFeature.H \
		$(TIGRPLUSPLUS)/TigrGffFeature.C
	g++ $(CFLAGS) -o $(OBJ)/TigrGffFeature.o -c \
		$(TIGRPLUSPLUS)/TigrGffFeature.C

baum-welch: \
		$(OBJ)/regex.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/baum-welch.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelch.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/Symbol.o
	g++ $(LDFLAGS) -o baum-welch \
		$(OBJ)/regex.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/baum-welch.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelch.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/Symbol.o \
		$(LIBS)

deterministic-trainer: \
		$(OBJ)/regex.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/deterministic-trainer.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelch.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/Symbol.o
	g++ $(LDFLAGS) -o deterministic-trainer \
		$(OBJ)/regex.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/deterministic-trainer.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelch.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/Symbol.o \
		$(LIBS)

model-combiner: \
		$(OBJ)/regex.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/model-combiner.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelch.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/Symbol.o
	g++ $(LDFLAGS) -o model-combiner \
		$(OBJ)/regex.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/model-combiner.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelch.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/Symbol.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/unveil.o:\
		$(UNVEIL)/unveil.C
	$(CC) $(CFLAGS) -o $(OBJ)/unveil.o -c \
		$(UNVEIL)/unveil.C
#--------------------------------------------------------
$(OBJ)/HMMGraph.o:\
		$(UNVEIL)/HMMGraph.C \
		$(UNVEIL)/HMMGraph.H
	$(CC) $(CFLAGS) -o $(OBJ)/HMMGraph.o -c \
		$(UNVEIL)/HMMGraph.C
#--------------------------------------------------------
$(OBJ)/FastForward.o:\
		$(UNVEIL)/FastForward.C \
		$(UNVEIL)/FastForward.H
	$(CC) $(CFLAGS) -o $(OBJ)/FastForward.o -c \
		$(UNVEIL)/FastForward.C
#---------------------------------------------------------
unveil: \
		$(OBJ)/regex.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/TigrGffReader.o \
		$(OBJ)/TigrGffExon.o \
		$(OBJ)/TigrGffFeature.o \
		$(OBJ)/TigrGffTranscript.o \
		$(OBJ)/TigrExceptions.o \
		$(OBJ)/TigrBitSet.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/unveil.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/Symbol.o
	$(CC) $(LDFLAGS) -o unveil \
		$(OBJ)/regex.o \
		$(OBJ)/TigrGffReader.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/TigrGffExon.o \
		$(OBJ)/TigrGffFeature.o \
		$(OBJ)/TigrGffTranscript.o \
		$(OBJ)/TigrExceptions.o \
		$(OBJ)/TigrBitSet.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/unveil.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/Symbol.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/unveil-eval-path.o:\
		$(UNVEIL)/unveil-eval-path.C
	$(CC) $(CFLAGS) -o $(OBJ)/unveil-eval-path.o -c \
		$(UNVEIL)/unveil-eval-path.C
#---------------------------------------------------------
unveil-eval-path: \
		$(OBJ)/TigrGffReader.o \
		$(OBJ)/TigrGffExon.o \
		$(OBJ)/TigrGffFeature.o \
		$(OBJ)/TigrGffTranscript.o \
		$(OBJ)/TigrExceptions.o \
		$(OBJ)/TigrBitSet.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/unveil-eval-path.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/Symbol.o
	$(CC) $(LDFLAGS) -o unveil-eval-path \
		$(OBJ)/TigrGffReader.o \
		$(OBJ)/TigrGffExon.o \
		$(OBJ)/TigrGffFeature.o \
		$(OBJ)/TigrGffTranscript.o \
		$(OBJ)/TigrExceptions.o \
		$(OBJ)/TigrBitSet.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/unveil-eval-path.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/Symbol.o \
		$(LIBS)


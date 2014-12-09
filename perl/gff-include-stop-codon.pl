#!/usr/bin/perl
#################################################################
# gff-include-stop-codon.pl
# bmajoros@tigr.org, Feb 2004
# This software is govered by the Artistic License.
#################################################################
use strict;
use SummaryStats;
use Translation;
use FastaWriter;
use FastaReader;
use FileHandle;
use GffTranscriptReader;
use Transcript;
$|=1;

use Getopt::Std;
our $opt_f;
getopts('f');

my $usage="$0 <in.gff> <*.fasta> <out.gff> [-f]
where -f=force final exon extension, even if no stop codon found";
die "$usage\n" unless @ARGV==3;
my ($gffFilename,$seqFilename,$outFilename)=@ARGV;

my %isStop; $isStop{"TAG"}=$isStop{"TGA"}=$isStop{"TAA"}=1;

# Read the sequences from the FASTA file
my %sequences;
my $fastaReader=new FastaReader($seqFilename);
while(1)
  {
    my ($defline,$seqRef)=$fastaReader->nextSequenceRef();
    last unless defined $defline;
    $defline=~/>\s*(\S+)/;
    my $id=$1;
    $sequences{$id}=$seqRef;
  }
$seqFilename=~/([^\/]+)$/;
my $seqBase=$1;

# Read the feature coordinates from the GFF file
my $gffReader=new GffTranscriptReader;
$gffReader->doNotSortTranscripts();
my $transcripts=$gffReader->loadGFF($gffFilename);

my $fastaWriter=new FastaWriter;
my $filehandle=new FileHandle(">$outFilename") ||
    die "Can't create $outFilename";
my $gff;
my $numTranscripts=@$transcripts;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $trans=$transcripts->[$i];
    my $substrate=$trans->{substrate};
    my $genomicSequence=$sequences{$substrate} || die $substrate;
    my $substrateLen=length $genomicSequence;
    my $transSeq=$trans->loadTranscriptSeq($genomicSequence);
    my $transcriptId=$trans->{transcriptId};
    my $numExons=@{$trans->{exons}};
    my $strand=$trans->{strand};
    my $begin=$trans->{begin};
    my $end=$trans->{end};
    my $stopLexeme=substr($transSeq,length($transSeq)-3,3);
    my $lastExon=$trans->getIthExon($numExons-1);
    my $lastExonType=$lastExon->getType();
    my $firstExon=$trans->getIthExon(0);
    my $firstExonType=$firstExon->getType();
    my $firstCodon=substr($transSeq,0,3);
    if(($firstExonType eq "initial-exon" || $firstExonType eq "single-exon")
       && $firstCodon ne "ATG")
      {
	print "WARNING! $firstCodon is not ATG...deleting transcript\n";
	next;
      }
    if(($lastExonType eq "final-exon" || $lastExonType eq "single-exon") 
       && (!$isStop{$stopLexeme} || $opt_f))
      {
	if($strand eq "+")
	  {if($lastExon->{end}+3<=$substrateLen) {$lastExon->{end}+=3}}
	else
	  {if($lastExon->{begin}-3>=0) {$lastExon->{begin}-=3}}
	undef $trans->{sequence};
	my $len=length $transSeq;
	undef $trans->getIthExon(0)->{sequence};
	$transSeq=$trans->loadTranscriptSeq($genomicSequence);
	my $newLen=length $transSeq;
	$stopLexeme=substr($transSeq,length($transSeq)-3,3);
	if(!$isStop{$stopLexeme} && !$opt_f)
	  {
	    print "WARNING! $stopLexeme is not a stop codon; deleting transcript $transcriptId\n(exon type=$lastExonType)\n";
	    undef $transcripts->[$i];
	    next;
	  }
      }
    my $newGff=$trans->toGff();
    $gff.=$newGff;
    undef $trans->{sequence};
  }
close($filehandle);

open(OUT,">$outFilename") || die "can't write to file: $outFilename\n";
print OUT "$gff\n";
close(OUT);




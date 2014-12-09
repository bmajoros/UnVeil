#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/perlib','.','/home/bmajoros/genomics/perl');
use FastaReader;
use FastaWriter;

my $usage="$0 <infile> <outfile>";
die "$usage\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

open(OUT,">$outfile") || die "can't create $outfile";
my $writer=new FastaWriter;

my $reader=new FastaReader($infile);
while(1)
  {
    my ($defline,$seq)=$reader->nextSequence();
    last unless defined $defline;
    scramble(\$seq);
    $writer->addToFasta($defline,$seq,\*OUT);
  }
close(OUT);

sub scramble
  {
    my ($seq)=@_;
    my $len=length $$seq;

    for(my $i=0 ; $i<$len ; ++$i)
      {
	my $target=int(rand($len));
	swap($seq,$i,$target);
      }
  }

sub swap
  {
    my ($seq,$x,$y)=@_;
    my $xChar=substr($$seq,$x,1);
    substr($$seq,$x,1)=substr($$seq,$y,1);
    substr($$seq,$y,1)=$xChar;
  }


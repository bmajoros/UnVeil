#!/usr/bin/perl
use strict;
use FastaReader;

die "$0 <*.fasta>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $reader=new FastaReader($infile);
while(1)
  {
    my ($def,$seq)=$reader->nextSequence();
    last unless defined $def;
    my $gc=($seq=~s/([GC])/$1/g);
    my $gcat=($seq=~s/([GCAT])/$1/g);
    my $content=int(100*$gc/$gcat+0.5);
    print "$content\n";
  }

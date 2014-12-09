#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;

my $usage="$0 \"pattern\" <*.fasta>";
die "$usage\n" unless @ARGV==2;
my ($pattern,$filename)=@ARGV;

my $writer=new FastaWriter();
my $reader=new FastaReader($filename);
while(1)
  {
    my ($defline,$sequence)=$reader->nextSequence();
    last unless defined $defline;
    next unless $defline=~/$pattern/;
    $writer->addToFasta($defline,$sequence,\*STDOUT);
  }



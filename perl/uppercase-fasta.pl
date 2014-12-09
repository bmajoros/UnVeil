#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/genomics/perl','/home/bmajoros/perlib');
use FastaReader;
use FastaWriter;
use TempFilename;

my $usage="$0 <filename>";
die "$usage\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $tmpFile=TempFilename::generate();
my $reader=new FastaReader($infile);
my $writer=new FastaWriter;
open(OUT,">$tmpFile");
while(1)
  {
    my ($defline,$sequence)=$reader->nextSequence();
    last unless defined $defline;
    $sequence="\U$sequence";
    if($defline=~/\n/) {chop $defline}
    if($sequence=~/\n/) {chop $sequence}
    $writer->addToFasta($defline,$sequence,\*OUT);
    #print OUT "$defline\n$sequence\n";
  }
close(OUT);

system("mv $tmpFile $infile");


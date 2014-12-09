#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/genomics/perl','/home/bmajoros/perlib');
use FastaWriter;

my $usage="cat seq.txt | $0 <identifier>";
die "$usage\n" unless @ARGV==1;
my ($id)=@ARGV;

my $seq;
while(<STDIN>)
  {
    $_=~s/\s+//g;
    $seq.=$_;
  }

my $writer=new FastaWriter;
$writer->addToFasta(">$id",$seq,\*STDOUT);


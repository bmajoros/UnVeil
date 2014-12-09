#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/perlib',
#	'/home/bmajoros/genomics/perl');
use FastaReader;

my $usage="$0 <*.fasta> <seq-id-or-dot(.)> <begin> <end>\n(coordinates are zero-based/space-based)";
die "$usage\n" unless @ARGV==4;
my ($fasta,$seqId,$begin,$end)=@ARGV;
my $length=$end-$begin;

my $reader=new FastaReader($fasta);
while(1)
  {
    my ($defline,$seq)=$reader->nextSequence();
    last unless defined($defline);
    if($seqId ne ".")
      {
	$defline=~/$\s*>\s*([^\s;]+)/ || die;
	my $id=$1;
	next unless $id eq $seqId;
      }
    my $subseq=substr($seq,$begin,$length);
    print "$subseq\n";
    exit unless $seqId eq ".";
  }



#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/perlib',
#	'/home/bmajoros/genomics/perl');
use FastaReader;

my $usage="$0 <*.fasta>";
die "$usage\n" unless @ARGV==1;
my ($filename)=@ARGV;

my ($numAs,$numTs,$numCs,$numGs);
my $reader=new FastaReader($filename);
my $numSeqs;
my $lenSum;
while(1)
  {
    my ($defline,$sequence)=$reader->nextSequence();
    last unless defined $defline;

    $numAs+=($sequence=~s/A/A/g);
    $numTs+=($sequence=~s/T/T/g);
    $numCs+=($sequence=~s/C/C/g);
    $numGs+=($sequence=~s/G/G/g);
    ++$numSeqs;
    $lenSum+=length($sequence);
  }

my $total=$numAs+$numTs+$numCs+$numGs;
my $percentA=int(10000*$numAs/$total)/100;
my $percentT=int(10000*$numTs/$total)/100;
my $percentC=int(10000*$numCs/$total)/100;
my $percentG=int(10000*$numGs/$total)/100;

my $pA=$numAs/$total;
my $pT=$numTs/$total;
my $pC=$numCs/$total;
my $pG=$numGs/$total;

my $H=-($pA*lg($pA) + 
	$pT*lg($pT) +
	$pC*lg($pC) +
	$pG*lg($pG));
my $Hmax=-lg(1/4.0);
my $Hnorm=$H/$Hmax;

my $percentGC=$percentG+$percentC;
my $percentAT=$percentA+$percentT;

my $meanLen=int($lenSum/$numSeqs+0.5);

print "$total bp\n";
print "$percentA% A, $percentT% T, $percentC% C,$percentG% G\n";
print "$percentGC% GC, $percentAT% AT\n";
print "H=$H H/Hmax=$Hnorm\n";
print "num seqs: $numSeqs    mean len=$meanLen\n";

sub lg
  {
    my ($x)=@_;
    return log($x)/log(2);
  }

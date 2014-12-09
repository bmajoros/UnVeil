#!/usr/bin/perl
use strict;
use GenbankParser;
use FastaWriter;

$0=~/([^\/]+$)/;
my $program=$1;
my $usage="$program <infile> <gff-outfile> <fasta-outfile>\n";
die $usage unless @ARGV==3;
my ($infile,$gffFile,$fastaFile)=@ARGV;

my $fastaWriter=new FastaWriter;
open(GFF,">$gffFile") || die "Can't create file: $gffFile\n";
open(FASTA,">$fastaFile") || die "Can't create file: $fastaFile\n";
my $transgrp=1;
my $parser=new GenbankParser($infile);
while(1)
  {
    my $entry=$parser->nextEntry();
    last unless $entry;

    my $defline=">$transgrp";
    my $organism=$entry->findUnique("ORGANISM") || 
      die "No ORGANISM field in Genbank entry!\n";
    my @words=split/\s+/,$organism;
    $organism="$words[0]_$words[1]";
    my $seq=$entry->getSubstrate();

    my $exons=$entry->getCDS();
    die "No CDS in Genbank entry!\n" unless $exons;
    my $numExons=@$exons;
    for(my $i=0 ; $i<$numExons ; ++$i)
      {
	my $exon=$exons->[$i];
	my ($begin,$end)=@$exon;
	my $strand=($begin<$end ? "+" : "-");
	if($strand eq "-") {($begin,$end)=($end,$begin)}
	print GFF "$transgrp\t$organism\texon\t$begin\t$end\t.\t$strand\t.\ttransgrp=$transgrp\n";
	$fastaWriter->addToFasta($defline,$seq,\*FASTA);
      }
    ++$transgrp;
  }
close(FASTA);
close(GFF);





package GffReader;
use strict;
use Exon;
use Feature;

######################################################################
#
# GffReader.pm 
# bmajoros@tigr.org 7/8/2002
#
# Returns a list of Features, sorted by their begin coordinates.
#
# Attributes:
#
# Methods:
#   $reader=new GffReader();
#   $featureArray=$reader->loadGFF($filename);
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
sub new
{
  my ($class)=@_;
  
  my $self={};
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#----------------------------------------------------------------
# $featureArray=$reader->loadGFF($filename);
#
sub loadGFF
  {
    my ($self,$gffFilename)=@_;
    my @features;
    if($gffFilename=~/\.gz$/)
      {open(GFF,"cat $gffFilename|gunzip|") || die $gffFilename}
    else
      {open(GFF,$gffFilename) || die $gffFilename}
    while(<GFF>)
      {
	next unless $_=~/\S+/;
	next if $_=~/^\s*\#/;
	my @fields=split/\s+/,$_;
	next unless @fields>7;
	my @additionalFields=splice(@fields,9);
	my $feature=new Feature($fields[0],$fields[1],$fields[2],
				$fields[3]-1,$fields[4],$fields[5],
				$fields[6],$fields[7],$fields[8],
				\@additionalFields);
	push @features,$feature;
      }
    close(GFF);
    @features=sort {$a->{begin} <=> $b->{begin}} @features;
    return \@features;
  }
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
#--------------------------------------------------------------------------

1;


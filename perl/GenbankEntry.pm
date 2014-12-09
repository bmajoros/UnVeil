package GenbankEntry;
use strict;

######################################################################
#
# GenbankEntry.pm bmajorostigr.org 11/16/2004
#
# Just an array of (key,value) pairs, where each value is either a
# string, a number, or a nested GenbankEntry.  NOTE that the keys are
# *not* assumed to be unique!
#
#
# Attributes:
#
# Methods:
#   $entry=new GenbankEntry();
#   $entry->addPair($key,$value);
#   my $n=$entry->numPairs();
#   my $pair=$entry->getIthPair($i); # returns pointer to [key,value]
#   my $pairs=$entry->findPairs($key); # returns array of pairs
#   my $value=$entry->findUnique($key); # returns one value element
#   $entry->print(*STDOUT);
#   my $cds=$entry->getCDS(); # returns array of [begin,end] pairs
#   my $seq=$entry->getSubstrate();
#
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $entry=new GenbankEntry();
sub new
{
  my ($class)=@_;

  my $self=
    {
     array=>[]
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $entry->addPair($key,$value);
sub addPair
  {
    my ($self,$key,$value)=@_;
    my $array=$self->{array};
    push @$array,[$key,$value];
  }
#---------------------------------------------------------------------
#   my $n=$entry->numPairs();
sub numPairs
  {
    my ($self)=@_;
    my $array=$self->{array};
    my $n=@$array;
    return $n;
  }
#---------------------------------------------------------------------
#   my $pair=$entry->getIthPair($i); # returns pointer to [key,value]
sub getIthPair
  {
    my ($self,$i)=@_;
    my $array=$self->{array};
    return $array->[$i];
  }
#---------------------------------------------------------------------
#   my $pairs=$entry->findPairs($key); # returns array of pairs
sub findPairs
  {
    my ($self,$key)=@_;
    my $array=$self->{array};
    my $n=@$array;
    my $pairs=[];
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $pair=$array->[$i];
	if($pair->[0] eq $key) {push @$pairs,$pair}
      }
    return $pairs;
  }
#---------------------------------------------------------------------
#   my $value=$entry->findUnique($key); # returns one value element
sub findUnique
  {
    my ($self,$key)=@_;
    my $array=$self->{array};
    my $n=@$array;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $pair=$array->[$i];
	if($pair->[0] eq $key) {return $pair->[1]}
      }
    return undef;
  }
#---------------------------------------------------------------------
sub print
  {
    my ($self,$handle,$level)=@_;
    my $array=$self->{array};
    my $n=@$array;
    my $pad=' 'x($level*3);
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $pair=$array->[$i];
	my ($key,$value)=@$pair;
	print "$pad$key => ";
	if(UNIVERSAL::isa($value,"GenbankEntry"))
	  {
	    print "\n$pad   \{\n";
	    $value->print($handle,$level+1);
	    print "$pad   }\n";
	  }
	else {print "$value\n"}
      }
  }
#---------------------------------------------------------------------
#   my $cds=$entry->getCDS(); # returns array of [begin,end] pairs
sub getCDS
  {
    my ($self)=@_;
    my $features=$self->findUnique("FEATURES");
    die "No FEATURES element found in Genbank entry!\n" unless $features;
    my $seq=$self->getSubstrate();
    my $len=length($seq);
    my $rhs=$features->findUnique("CDS");
    die "No CDS element found in Genbank FEATURES clause!\n"
      unless $rhs;
    my $array=[];
    if($rhs=~/^\s*(complement\()?([A-Za-z0-9]+:)?(\d+)\.\.(\d+)/)
      {
	my ($complement,$junk,$begin,$end)=($1,$2,$3,$4);
	if($junk)### TOTALLY BOGUS...FOR BURGE'S EST FILE
	  {
	    push @$array,[1,$len];
	    return $array;
	  }
	if($complement) {($begin,$end)=($end,$begin)}
	push @$array,[$begin,$end];
      }
    elsif($rhs=~/^\s*(complement\()?join\(([A-Za-z0-9]+:)?([^)]+)\)/)
      {
	my ($complement,$junk,$coords)=($1,$2,$3);
	if($junk)### TOTALLY BOGUS...FOR BURGE'S EST FILE
	  {
	    push @$array,[1,$len];
	    return $array;
	  }
	$coords=~s/\s+//g;
	my @exons=split/,/,$coords;
	foreach my $exon (@exons)
	  {
	    $exon=~/(\d+)\.\.(\d+)/ || 
	      die "Can't parse Genbank exon: $exon\n";
	    my ($begin,$end)=($1,$2);
	    if($complement) {($begin,$end)=($end,$begin)}
	    push @$array,[$begin,$end];
	  }
      }
    else
      {die "Can't parse CDS clause in Genbank entry: $rhs\n"}
    return $array;
  }
#---------------------------------------------------------------------
#   my $seq=$entry->getSubstrate();
sub getSubstrate
  {
    my ($self)=@_;
    my $seq=$self->findUnique("ORIGIN");
    $seq="\U$seq";
    return $seq;
  }
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;


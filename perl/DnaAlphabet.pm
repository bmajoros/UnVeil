package DnaAlphabet;
use strict;
use Alphabet;
use vars qw(@ISA);

@ISA=qw(Alphabet);

######################################################################
#
# DnaAlphabet.pm bmajoros@tigr.org 10/10/2003
#
# 
# 
#
# Attributes:
#
# Methods:
#   $dnaAlphabet=new DnaAlphabet();
#
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
sub new
{
  my ($class)=@_;

  my $self=new Alphabet("ACGNT");
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;


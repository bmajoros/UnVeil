package GenbankParser;
use strict;
use FileHandle;
use GenbankEntry;

######################################################################
#
# GenbankParser.pm bmajorostigr.org 11/15/2004
#
# 
# 
#
# Attributes:
#
# Methods:
#   $parser=new GenbankParser($filename);
#   $entry=$parser->nextEntry(); # returns a GenbankEntry object
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $parser=new GenbankParser($filename);
sub new
{
  my ($class,$filename)=@_;

  my $filehandle=new FileHandle($filename) || die "can't open $filename";
  my $self=
    {
     file=>$filehandle
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $entry=$parser->nextEntry();
sub nextEntry
  {
    my ($self)=@_;
    my $file=$self->{file};
    if(eof($file)) {return undef}
    my $entry=new GenbankEntry;
    while(<$file>)
      {
	# Each line is formatted into a 12-byte label field
	# followed by a varible-length part
      REDO:
	my $key=substr($_,0,12);
	$key=~s/\s+//g;
	my $value=substr($_,12,-1);
	$value=~s/^\s+//g;
	$value=~s/\s+$//g;
	if($key eq "FEATURES")
	  {
	    $value=new GenbankEntry;
	    $entry->addPair($key,$value);
	    while(<$file>)
	      {
		# Now the label field is 21-bytes		
		last unless $_=~/^\s/;
		my $label=substr($_,0,21);
		my $rhs=substr($_,21,-1);
		$rhs=~s/^\s+//g;
		$rhs=~s/\s+$//g;
		$label=~s/\s+//g;
		if(length($label)>0)
		  {$value->addPair($label,$rhs)}
		else
		  {
		    my $pair=$value->getIthPair($value->numPairs()-1);
		    $pair->[1].=" $rhs";
		  }
	      }
	    goto REDO;
	  }
	elsif($key eq "ORIGIN")
	  {
	    $value="";
	    while(<$file>)
	      {
		last if(/\/\//);
		$_=~s/[\d\s]+//g;
		$value.=$_;
	      }
	  }
	if(length($key)>0)
	  {$entry->addPair($key,$value)}
	else
	  {
	    my $pair=$entry->getIthPair($entry->numPairs()-1);
	    $pair->[1].=" $value";
	  }
	last if($key eq "ORIGIN");
      }
    return $entry;
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;


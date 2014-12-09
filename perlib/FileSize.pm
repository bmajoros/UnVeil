package FileSize;
use strict;

######################################################################
#
# FileSize.pm bmajoros 2/27/2001
#
# Functions:
#   $fileSize=FileSize::fileSize($filename);
#
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
sub fileSize
  {
    my ($filename)=@_;
    my @stats=stat $filename;
    return $stats[7];
  }
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;


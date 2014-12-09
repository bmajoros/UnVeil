package SummaryStats;
use strict;
use Carp;

######################################################################
#
# SummaryStats.pm bmajoros 2/19/2001
#
# 
# 
# ($mean,$stddev,$min,$max)=SummaryStats::summaryStats(\@array);
# ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats(\@array);
# $r=correlation(\@array1,\@array2);
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
# ($mean,$stddev,$min,$max)=summaryStats(\@array);
sub summaryStats
  {
    my ($array)=@_;
    my $n=@$array;
    my ($minX,$maxX);
    my $sumX=0;
    my $sumXX=0;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $x=$array->[$i];
	$sumX+=$x;
	$sumXX+=($x*$x);
	$minX=$maxX=$x if $i==0;
	$minX=$x if $x<$minX;
	$maxX=$x if $x>$maxX;
      }
    my $meanX=$sumX/$n;
    #confess "n\=$n in summaryStats()\n" unless $n>1;
    my $varX=$n>1 ? ($sumXX-$sumX*$sumX/$n)/($n-1) : undef;
    my $stddevX=sqrt($varX);
    return ($meanX,$stddevX,$minX,$maxX);
  }
#---------------------------------------------------------------------
# $r=correlation(\@array1,\@array2);
sub correlation
  {
    my ($Xs,$Ys)=@_;

    my $sumX=0.0;
    my $sumY=0.0;
    my $sumXY=0.0;
    my $sumXX=0.0;
    my $sumYY=0.0;
    my $n=@$Xs;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $x=$Xs->[$i];
	my $y=$Ys->[$i];

	$sumX+=$x;
	$sumY+=$y;
	$sumXY+=($x*$y);
	$sumXX+=($x*$x);
	$sumYY+=($y*$y);
      }
    
    my $r=($sumXY-$sumX*$sumY/$n)/
      sqrt(($sumXX-$sumX*$sumX/$n)*($sumYY-$sumY*$sumY/$n));
    return $r;
  }
#---------------------------------------------------------------------
# ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats(\@array);
sub roundedSummaryStats
  {
    my ($array)=@_;
    my ($mean,$stddev,$min,$max)=SummaryStats::summaryStats($array);
    $mean=int(10*$mean+0.5)/10;
    $stddev=int(10*$stddev+0.5)/10;
    $min=int(10*$min+0.5)/10;
    $max=int(10*$max+0.5)/10;
    return ($mean,$stddev,$min,$max);
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------





#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;


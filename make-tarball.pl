#!/usr/bin/perl
use strict;

my $version="1.0";

my $tarfile="unveil-src-$version.tar";
if(-e $tarfile) {unlink($tarfile)}
my $gz="${tarfile}.gz";
if(-e $gz) {unlink($gz)}
my $command="tar cf $tarfile README unveil model-combiner deterministic-trainer baum-welch *.[HC] makefile *.pl tigr++/*.[HC] doc/training.txt perl/*.pl perl/*.pm perlib/*.pm metamodel.hmms submodels.txt";
system("$command");
`gzip $tarfile`;

print "packaged into $tarfile.gz\n";


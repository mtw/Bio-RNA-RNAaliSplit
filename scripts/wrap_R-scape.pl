#!/usr/bin/env perl
# Last changed Time-stamp: <2017-07-05 18:27:42 mtw>
# -*-CPerl-*-
#
# A wrapper for R-scape
#
# usage: wrap_R-scape.pl -a myaln.stk --statistics RAFS
#

use version; our $VERSION = qv('0.05.2');
use strict;
use warnings;
use Bio::RNA::RNAaliSplit::WrapRscape;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Data::Dumper;
use Pod::Usage;
use Path::Class;
use Carp;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my $show_version = 0;
my ($id,$stkfile,$stat);

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("a|aln=s"        => \$stkfile,
					   "s|statistic=s" => \$stat,
                                           "version"     => sub{$show_version = 1},
                                           "man"         => sub{pod2usage(-verbose => 2)},
                                           "help|h"      => sub{pod2usage(1)}
                                           );

if ($show_version == 1){
  print "wrap_R-scape.pl $VERSION\n";
  exit(0);
}

unless (-f $stkfile){
  warn "Could not find input file provided via -a|-aln option";
  pod2usage(-verbose => 0);
}

my $r = Bio::RNA::RNAaliSplit::WrapRscape->new(ifile => $stkfile,
					       statistic => $stat,
					       odir => ['test'],
					      );
print Dumper($r);

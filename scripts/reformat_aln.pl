#!/usr/bin/env perl
# Last changed Time-stamp: <2017-07-03 16:37:37 mtw>
# -*-CPerl-*-
#
# inter-convert alignment formats

use strict;
use warnings;
use Bio::AlignIO;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;

my $infile_aln = undef;
my $infmt = undef;
my $outfmt = undef;

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("a|aln=s"    => \$infile_aln,
					   "i|infmt=s"  => \$infmt,
					   "o|otrfmt=s" => \$outfmt,
                                           "man"        => sub{pod2usage(-verbose => 2)},
                                           "help|h"     => sub{pod2usage(1)}
					  );

unless (-f $infile_aln){
  warn "Could not find input alignment file provided via -a|--aln option";
  pod2usage(-verbose => 0);
}

my $in = Bio::AlignIO->new(-file => "$infile_aln" ,
			   -format => $infmt);

my $out = Bio::AlignIO->new(-fh   => \*STDOUT ,
			 -format => $outfmt);

while ( my $aln = $in->next_aln ) {
    $out->write_aln($aln);
}


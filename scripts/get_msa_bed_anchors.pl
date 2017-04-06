#!/usr/bin/env perl
# Last changed Time-stamp: <2017-04-06 11:28:29 mtw>
# -*-CPerl-*-
#
# Create BED6 anchors for columns in a MSA
#
# N.B.: Input aln file MUST NOT contain sequence limits !!!
# (e.g. ZIKV_AFR3/1890-2340) but only the display_id, e.g. ZIKV_AFR3
#
# usage: get_msa_bed_anchors.pl -a myaln.aln -c 20,133,176
#

use version; our $VERSION = qv('0.05_01');
use strict;
use warnings;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::Location::Simple;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Data::Dumper;
use Pod::Usage;
use Path::Class;
use Carp;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my $format = "ClustalW";
my $show_version = 0;
my ($id,$alifile,$columns,$anchor_name);
my $anchor_nr = 1;

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("aln|a=s"     => \$alifile,
					   "columns|c=s" => \$columns,
					   "version"     => sub{$show_version = 1},
                                           "man"         => sub{pod2usage(-verbose => 2)},
                                           "help|h"      => sub{pod2usage(1)}
                                           );

if ($show_version == 1){
  print "get_msa_bed_anchors.pl $VERSION\n";
  exit(0);
}

unless (-f $alifile){
  warn "Could not find input file provided via --aln|-a option";
  pod2usage(-verbose => 0);
}

my @cols = split /,/, $columns;

croak "[ERROR] must provide a comman-separated list of --columns|-c option"
  unless (scalar(@cols) >= 1);
foreach (@cols){
  croak "[ERROR] argument $_ of --columns|-c is not a number"
    unless (/^\d+$/);
}

my $input_AlignIO = Bio::AlignIO->new(-file => $alifile,
				      -format => 'ClustalW'
				     );

my $input_aln = $input_AlignIO->next_aln;

open my $out, ">", "anchors.bed" or die $!;

foreach my $c (@cols){
  $anchor_name = join "-", ("a",$anchor_nr,$c);
  foreach my $seq ($input_aln->each_seq) {
    my $loc = $seq->location_from_column($c);
    my $id = $seq->display_id;
    croak "[ERROR] column $c not defined for sequence $id, please choose larger column number"
      unless (defined($loc));
    croak "[ERROR] invald position: $id column $c"
      unless ($loc->valid_Location);
    croak "[ISSUE] column $c is a gap (type IN-BETWEEN) in sequence $id"
      if ($loc->location_type() eq "IN-BETWEEN");
    # print Dumper ($loc);
    my $start = $loc->start;
    #    print "$id $start\n";
    my $FeatureIntervalN = join "\t", ($id, $start-1, $start, $anchor_name);
    print $out $FeatureIntervalN."\n";
  }
  $anchor_nr++;
  #print "------------------------\n";
}

close $out;
#print Dumper($anchors);

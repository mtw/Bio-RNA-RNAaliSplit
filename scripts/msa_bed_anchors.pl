#!/usr/bin/env perl
# Last changed Time-stamp: <2019-04-05 22:58:09 mtw>
# -*-CPerl-*-
#
# Create BED6 anchors for columns in a multiple sequence alignment
#
# N.B.: Input aln file MUST NOT contain sequence limits !!!
# (e.g. ZIKV_AFR3/1890-2340) but only the display_id, e.g. ZIKV_AFR3
#
# usage: msa_bed_anchors.pl -a myaln.aln -c 20,133,176
#

use version; our $VERSION = qv('0.10');
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
my @regions=();

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("aln|a=s"     => \$alifile,
					   "columns|c=s" => \$columns,
					   "version"     => sub{$show_version = 1},
					   "man"         => sub{pod2usage(-verbose => 2)},
					   "help|h"      => sub{pod2usage(1)}
					   );

if ($show_version == 1){
  print "msa_bed_anchors.pl $VERSION\n";
  exit(0);
}

unless (-f $alifile){
  warn "Could not find input file provided via --aln|-a option";
  pod2usage(-verbose => 0);
}

@regions = split /,/, $columns;

croak "[ERROR] must provide a comman-separated list to --columns|-c option"
  unless (scalar(@regions) >= 1);
foreach (@regions){
  print ">>$_<<\n";
  unless (/^\d+$/){
    unless (/^\d+\-\d+$/){
      croak "[ERROR] argument $_ of --columns|-c is not a (int-int) region";
    }
    else{
      next;
    }
    croak "[ERROR] argument $_ of --columns|-c is not a int number"
  }
}

print Dumper(\@regions);

my $input_AlignIO = Bio::AlignIO->new(-file => $alifile,
				      -format => $format,
				     );

my $input_aln = $input_AlignIO->next_aln;

open my $out, ">", "anchors.bed" or die $!;

foreach my $c (@regions){
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
}

close $out;


__END__

=head1 NAME

msa_bed_anchors.pl - Create BED anchors for multiple sequence alignments

=head1 SYNOPSIS

msa_bed_anchors.pl [--aln|-a I<FILE>] [--columns|-c I<LIST>] [options]

=head1 DESCRIPTION

This tool creates BED anchors from columns or regions in a multiple
sequence alignment (MSA) that can be used to produce a
anchor-constrained MSA. A typical use case is refinement of Clustal or
MAFFT alignments with anchored mlocarna alignments.

Anchor columns or regions relative to the input MSA are transformed
into BED-type intervals for each sequence of the MSA. An BED output
file 'anchors.bed' that can directly be passed to the mlocarna
--anchor-constraints option is written to the current working
directory.

=head1 OPTIONS

=over

=item B<--aln|-a>

A multiple sequence alignment in ClustalW format.
IMPORTANT: The sequence identifiers MUST NOT contain sequence limits,
e.g. ZIKV_AFR3/1890-2340 but only the display_id, e.g. ZIKV_AFR3 .

=item B<--columns|-c>

Comma-separated list of columns (int) or regions (int-int) in the
input MSA that correspond to anchor constraints for the new
alignment. A single column will result in an anchor constraint of
length 1, whereas a region will give an anchor constraint of respective
length in the output BED file.

=item B<--version>

Show msa_bed_anchors.pl version and exit

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt> and
E<lt>michael.wolfinger@univie.ac.atE<gt>

=cut

#!/usr/bin/env perl
# Last changed Time-stamp: <2016-12-31 23:53:23 mtw>
# -*-CPerl-*-

use AlignSplit;
use Data::Dumper;
#use Bio::AlignIO;

my $alnfile = "./result.aln";
my $stkfile = "./result.stk";

my $format = "ClustalW";

my $AlignSplitObject = AlignSplit->new(file => $alnfile,
				       format => $format,
    );

#print Dumper($AlignSplitObject);

my $aln = $AlignSplitObject->alignment->next_aln();

# Extract sequences and check values for the alignment column $pos
foreach $seq ($aln->each_seq) {
 print $seq->seq()."\n";
}


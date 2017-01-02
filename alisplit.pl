#!/usr/bin/env perl
# Last changed Time-stamp: <2017-01-02 13:46:31 mtw>
# -*-CPerl-*-

use AlignSplit;
use Data::Dumper;
#use Bio::AlignIO;
use Bio::Align::DNAStatistics;
use Bio::Tree::DistanceFactory;
use Bio::TreeIO;

my $alnfile = "./result.aln";
my $stkfile = "./result.stk";

my $format = "ClustalW";

my $AlignSplitObject = AlignSplit->new(file => $alnfile,
				       format => $format,
    );

#print Dumper($AlignSplitObject);

my $aln = $AlignSplitObject->alignment->next_aln();

my $dfactory = Bio::Tree::DistanceFactory->new(-method => 'NJ');
my $stats = Bio::Align::DNAStatistics->new;
my $treeout = Bio::TreeIO->new(-format => 'newick');

print Dumper($dfactory);
print Dumper($stats);

my $jcmatrix = $stats->distance(-align => $aln, 
				-method => 'F81');

print Dumper(jcmatrix);
# Extract sequences and check values for the alignment column $pos
#foreach $seq ($aln->each_seq) {
# print $seq->seq()."\n";
#}

#### TODO read in distance matrix produced by another program, such as AnalyseSequences

#my $parser = Bio::Matrix::IO->new(-format => 'phylip',
#                                    -file   => 'filename.dist');
#my $mat  = $parser->next_matrix;
#my $tree = $dfactory->make_tree($mat);
#$treeout->write_tree($tree);

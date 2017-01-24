#!/usr/bin/env perl
# Last changed Time-stamp: <2017-01-24 17:40:37 mtw>
# -*-CPerl-*-
#
# usage: alisplit.pl -a myfile.aln
#
# NB: Display ID handling in Bio::AlignIO is broken for Stockholm
# format. Use ClustalW format instead !!!

use strict;
use warnings;
use AlignSplit;
use WrapRNAz;
use WrapAnalyseDists;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Data::Dumper;
use Pod::Usage;
use Carp;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my $format = "ClustalW";
my $method = "dHn"; # SCI | dHn | dHx
my @nseqs=();
my ($i,$j,$dim,$alnfile, $Dfile);
my @pw_alns = ();
my @D = ();
my $check = 1;

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("aln|a=s"    => \$alnfile,
                                           "method|m=s" => \$method,
                                           "man"        => sub{pod2usage(-verbose => 2)},
                                           "help|h"     => sub{pod2usage(1)}
                                           );

unless (-f $alnfile){
  warn "Could not find input file provided via --aln|-a option";
  pod2usage(-verbose => 0);
}


my $AlignSplitObject = AlignSplit->new(ifile => $alnfile,
				       format => $format,
				       dump => 1);
#print Dumper($AlignSplitObject);
#die;
$dim = $AlignSplitObject->next_aln->num_sequences;

# extract all pairwise alignments
print STDERR "Extracting pairwise alignments ...\n";
for ($i=1;$i<$dim;$i++){
  for($j=$i+1;$j<=$dim;$j++){
    my $token = join "_", "pw",$i,$j;
    my ($sa_clustal,$sa_stockholm) = $AlignSplitObject->dump_subalignment("pairwise", $token, [$i,$j]);
    push @pw_alns, $sa_clustal->stringify;
  }
}

# initialize distance matrix
for($i=0;$i<$dim;$i++){
  for ($j=0;$j<$dim;$j++){
    $D[$dim*$i+$j] = 0.;
  }
}

# build distance matrix based on pairwise alignments
print STDERR "Constructing distance matrix based on pairwise alignments ...\n";
foreach my $ali (@pw_alns){
  #print "processing $ali ...\n";
  my $pw_aso = AlignSplit->new(ifile => $ali,
			       format => "ClustalW");
  #print Dumper($pw_aso);
  my ($i,$j) = sort split /_/, $pw_aso->infilebasename;

  if ($method eq "SCI"){
    # distance = -log( [normalized] pairwise SCI)
    my ($sci,$dsci);
    if($pw_aso->sci > 1){$sci = 1}
    elsif ($pw_aso->sci == 0.){$sci = 0.000001}
    else { $sci = $pw_aso->sci; }
    $dsci = -1*log($sci)/log(10);
    $D[$dim*($i-1)+($j-1)] =  $D[$dim*($j-1)+($i-1)] = $dsci;
    #  print "$i -> $j : $sci\t$dsci\n";
  }
  elsif ($method eq "dHn") {
    # Hamming dist with gaps replaced by Ns
    $D[$dim*($i-1)+($j-1)] =  $D[$dim*($j-1)+($i-1)] = $pw_aso->hammingdistN;
  }
  elsif ($method eq "dHx") {
    # TODO compute $pw_aso->hammingdistX in AlignSplit.pm
    # Hamming dist with all gap columns removed
    $D[$dim*($i-1)+($j-1)] =  $D[$dim*($j-1)+($i-1)] = $pw_aso->hammingdistX;
  }
  else {croak "Method $method not available ..please use SCI|dHn|dHx"}
}



# write matrix to file
print STDERR "Writing distance matrix to file ...\n";
if ($method eq "SCI"){ $Dfile = dump_matrix(\@D,$dim,1,1,"S")}
elsif ($method eq "dHn"){$Dfile = dump_matrix(\@D,$dim,1,1,"dHn")}
elsif ($method eq "dHx"){$Dfile = dump_matrix(\@D,$dim,1,1,"dHx")}
else { croak "Method $method not available ..please use SCI|dHn|dHx"}

# check triangle inequality
if ($check == 1){
  print STDERR "Checking triangle inequality ...\n";
  my $result = check_triangle($dim,\@D);
}

# compute Neighbor Joining tree and do split decomposition
print STDERR "Perform Split Decomposition ...\n";
my $sd = WrapAnalyseDists->new(ifile => $Dfile,
			       odir => $AlignSplitObject->odir,
			       basename => $AlignSplitObject->infilebasename);
print "Identified ".$sd->count." splits\n";

# run RNAz for the input alignment
my $rnaz = WrapRNAz->new(ifile => $alnfile,
			 odir => $AlignSplitObject->odir);
print join "\t", "#SVM prob","SCI","Z-score","hit","sequences","alignment\n";
print join "\t", $rnaz->P,$rnaz->sci,$rnaz->z,"O",$dim,$alnfile."\n";
print "-----------------------------------------\n";

# extract split sets and run RNAz on each of them
my $splitnr=1;
while (my $sets = $sd->pop()){
  my ($rnazo1,$rnazo2,$have_rnazo1,$have_rnazo2,$ao1,$ao2,$hint) = (0)x7;
  my ($sa1_c,$sa1_s,$sa2_c,$sa2_s); # subalignments in Clustal and Stockholm
  my $set1 = $$sets{S1};
  my $set2 = $$sets{S2};
  my $token = "split".$splitnr;
#  print "set1: @$set1\n";
#  print "set2: @$set2\n";
  ($sa1_c,$sa1_s) = $AlignSplitObject->dump_subalignment("splits", $token.".set1", $set1);
  ($sa2_c,$sa2_s) = $AlignSplitObject->dump_subalignment("splits", $token.".set2", $set2);
  if( scalar(@$set1) > 1){
    $rnazo1 = WrapRNAz->new(ifile => $sa1_c,
			    odir => $AlignSplitObject->odir);
    $have_rnazo1 = 1;
    ($rnazo1->P > $rnaz->P) ? ($hint = "*") : ($hint = " ");
    $ao1 = WrapRNAalifold->new(ifile => $sa1_c,
			       odir => $AlignSplitObject->odir);
#    print Dumper($ao1);
    print join "\t",$rnazo1->P,$rnazo1->sci,$rnazo1->z,$hint,scalar(@$set1),$sa1_c."\n";
  }
  if( scalar(@$set2) > 1){
    $rnazo2 = WrapRNAz->new(ifile => $sa2_c,
			    odir => $AlignSplitObject->odir);
    $have_rnazo2 = 1;
    ($rnazo2->P > $rnaz->P) ? ($hint = "*") : ($hint = " ");
    $ao2 = WrapRNAalifold->new(ifile => $sa2_c,
			       odir => $AlignSplitObject->odir);
#    print Dumper($ao2);
    print join "\t",$rnazo2->P,$rnazo2->sci,$rnazo2->z,$hint,scalar(@$set2),$sa2_c."\n";
  }
  $splitnr++;
}


###############
# subroutines #
###############


sub dump_matrix {
  my ($M, $d, $ad, $pd, $what) = @_;
  my ($i,$j,$info);
  my $ad_mx = "ld.mx"; # AnalyseDists lower diagoinal distance matrix
  my $pd_mx = "pd.dst"; # Phylip distance matrix

  if ($what eq "S"){$info="> S (SCI distance)"}
  elsif($what eq "dHn"){$info="> H (Hamming distance with gap Ns)"}
  elsif($what eq "dHx"){$info="> H (Hamming distance with gaps removed)"}
  else{$info="> H (unknown)"}

  if (defined $ad){ # print lower triangle matrix input for AanalyseDists
    open my $matrix, ">", $ad_mx or die $!;
    print $matrix $info,"\n";
    print $matrix "> X $d\n";
    for ($i=1;$i<$d;$i++){
      for($j=0;$j<$i;$j++){
	printf $matrix "%6.4f ", @$M[$d*$i+$j];
      }
      print $matrix "\n";
    }
    close $matrix;
  }
  else{ # print entire matrix
    for ($i=0;$i<$d;$i++){
      for($j=0;$j<$d;$j++){
	printf "%6.4f ", @$M[$d*$i+$j];
      }
      print "\n";
    }
  }
  if (defined $pd){ # print phylip distance matrix
    open my $matrix, ">", $pd_mx  or die $!;
    print $matrix "$d\n";
    for ($i=0;$i<$d;$i++){
      print $matrix eval($i+1)." ";
	for($j=0;$j<$d;$j++){
	  printf $matrix "%6.4f ", @$M[$d*$i+$j];
	}
      print $matrix "\n";
    }
    close $matrix;
  }
  return $ad_mx;
}

sub check_triangle {
  my ($d,$dref) = @_;
  my $count = 0;
  my @M = @$dref;
  my ($i,$j,$k,$d_ij,$d_jk,$d_ik);
  for ($i=0;$i<$d;$i++){
    for($j=0;$j<$d;$j++){
      $d_ij = $M[$d*$i+$j];
      for($k=0;$k<$d;$k++){
	$d_jk = $M[$d*$j+$k];
	$d_ik = $M[$d*$i+$k];
	$count++;
	croak "ERROR triangle inequation not satisfied i:$i j:$j k:$k"
	  unless ($d_ij + $d_jk >= $d_ik);
      }

    }
  }
  print STDERR "Checked $count triangles ...\n";
}

__END__

=head1 NAME

alisplit.pl - Split and decompose multiple sequence alignments

=head1 SYNOPSIS

alisplit.pl [--aln|-a I<FILE>] [--method|-m I<OPTION>] [options]

=head1 DESCRIPTION

This tool splits multiple sequence alignments horizontally, thereby
extracting sets of sequences that group together according to a
decision value. The most natural decision value is the RNAz SVM
RNA-classs probability.

A tree is computed from pairwise distances of sequences in the input
alignment and subsets of the alignment are derived by performing a
split decomposition of the matrix of pairwise distances. These
subsets/subalignments are then evaluated according to the same
decision value and a decision is made whether a subalignment performs
better than the original alignment.This can be used to discriminate
sequences that to not 'fit' in the input alignment.

=head1 OPTIONS

=over

=item B<--aln|-a>

A multiple sequence alignment in ClustalW format

=item B<--method|-m>

Method to compute pairwise ditances. Available options are 'dHn',
'dHx' and 'SCI'. The first and second computes pairwise proximity as
Hamming distance on the s=level of sequences. 'dHn' replaces gaps with
'N', whereas 'dHx' removes all gap columns (is not yet
implemented). 'SCI' computes the distance as 1-log(SCI), based on a
truncated strucure conservation index of two sequences. The latter,
however, is not a metric and therefore often results in negative
branch lengths in Neighbor Joining trees. Use with caution. [default:
'dHn']

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut
(END)

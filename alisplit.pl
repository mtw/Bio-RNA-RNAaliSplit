#!/usr/bin/env perl
# Last changed Time-stamp: <2017-01-31 12:57:49 mtw>
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
use Path::Class;
use Carp;
use RNA;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my $format = "ClustalW";
my $method = "dHn"; # SCI | dHn | dHx | dBp
my $outdir = "as";
my @nseqs=();
my ($dim,$alnfile);

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("aln|a=s"    => \$alnfile,
                                           "method|m=s" => \$method,
					   "out|o=s"    => \$outdir,
                                           "man"        => sub{pod2usage(-verbose => 2)},
                                           "help|h"     => sub{pod2usage(1)}
                                           );

unless (-f $alnfile){
  warn "Could not find input file provided via --aln|-a option";
  pod2usage(-verbose => 0);
}


my $AlignSplitObject = AlignSplit->new(ifile => $alnfile,
				       format => $format,
				       odirn => $outdir,
				       dump => 1);
#print Dumper($AlignSplitObject);
print Dumper(${$AlignSplitObject->next_aln}{_order});
#die;
$dim = $AlignSplitObject->next_aln->num_sequences;

my $dmfile = make_distance_matrix($AlignSplitObject,$method);

# compute Neighbor Joining tree and do split decomposition
print STDERR "Perform Split Decomposition ...\n";
my $sd = WrapAnalyseDists->new(ifile => $dmfile,
			       odir => $AlignSplitObject->odir,
			       basename => $AlignSplitObject->infilebasename);
print STDERR "Identified ".$sd->count." splits\n";

# run RNAalifold for the input alignment
my $alifold = WrapRNAalifold->new(ifile => $alnfile,
				  odir => $AlignSplitObject->odir);
# run RNAz for the input alignment
my $rnaz = WrapRNAz->new(ifile => $alnfile,
			 odir => $AlignSplitObject->odir);
print join "\t", "#RNAz SVM prob","hit","z-score","SCI RNAz","SCI aifold","sequences","alignment\n";
print join "\t", $rnaz->P,"0",$rnaz->z,$rnaz->sci,$alifold->sci,$dim,$alnfile."\n";
print "-------------------------------------------------------------\n";

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
    if($rnazo1->P > 1.1*$rnaz->P){$hint = "**"};
    if($rnazo1->P > 1.2*$rnaz->P){$hint = "***"};
    if($rnazo1->P > 1.3*$rnaz->P){$hint = "****"};
    $ao1 = WrapRNAalifold->new(ifile => $sa1_c,
			       odir => $AlignSplitObject->odir);
#    print Dumper($ao1);
    print join "\t",$rnazo1->P,$hint,$rnazo1->z,$rnazo1->sci,$ao1->sci,scalar(@$set1),$sa1_c."\n";
  }
  if( scalar(@$set2) > 1){
    $rnazo2 = WrapRNAz->new(ifile => $sa2_c,
			    odir => $AlignSplitObject->odir);
    $have_rnazo2 = 1;
    ($rnazo2->P > $rnaz->P) ? ($hint = "*") : ($hint = " ");
    if($rnazo2->P > 1.1*$rnaz->P){$hint = "**"};
    if($rnazo2->P > 1.2*$rnaz->P){$hint = "***"};
    if($rnazo2->P > 1.3*$rnaz->P){$hint = "****"};
    $ao2 = WrapRNAalifold->new(ifile => $sa2_c,
			       odir => $AlignSplitObject->odir);
#    print Dumper($ao2);
    print join "\t",$rnazo2->P,$hint,$rnazo2->z,$rnazo2->sci,$ao2->sci,scalar(@$set2),$sa2_c."\n";
  }
  $splitnr++;
}


###############
# subroutines #
###############

sub make_distance_matrix {
  my ($ASO,$m) = @_;
  my ($i,$j,$Dfile);
  my @pw_alns = ();
  my @D = ();
  my $check = 1;
  my $dim = $ASO->next_aln->num_sequences;

  # extract all pairwise alignments
  print STDERR "Extracting pairwise alignments ...\n";
  for ($i=1;$i<$dim;$i++){
    for($j=$i+1;$j<=$dim;$j++){
      my $token = join "_", "pw",$i,$j;
      my ($sa_clustal,$sa_stockholm) = $ASO->dump_subalignment("pairwise", $token, [$i,$j]);
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
    my $pw_aso = AlignSplit->new(ifile => $ali,
				 format => "ClustalW",
				 odirn => $outdir);
    my ($i,$j) = sort split /_/, $pw_aso->infilebasename;

    if ($m eq "SCI"){
      # distance = -log( [normalized] pairwise SCI)
      my ($sci,$dsci);
      if($pw_aso->sci > 1){$sci = 1}
      elsif ($pw_aso->sci == 0.){$sci = 0.000001}
      else { $sci = $pw_aso->sci; }
      $dsci = -1*log($sci)/log(10);
      $D[$dim*($i-1)+($j-1)] =  $D[$dim*($j-1)+($i-1)] = $dsci;
      #  print "$i -> $j : $sci\t$dsci\n";
    }
    elsif ($m eq "dHn") {
      # Hamming dist with gaps replaced by Ns
      $D[$dim*($i-1)+($j-1)] =  $D[$dim*($j-1)+($i-1)] = $pw_aso->hammingdistN;
    }
    elsif ($m eq "dHx") {
      # TODO compute $pw_aso->hammingdistX in AlignSplit.pm
      # Hamming dist with all gap columns removed
      $D[$dim*($i-1)+($j-1)] =  $D[$dim*($j-1)+($i-1)] = $pw_aso->hammingdistX;
    }
    elsif ($m eq "dBp") {
      my $so1 = $pw_aso->next_aln->get_seq_by_pos(1);
      my $so2 = $pw_aso->next_aln->get_seq_by_pos(2);
      my $seq1 = $so1->seq;
      my $seq2 = $so2->seq;
      my ($ss1,$mfe1) = RNA::fold($seq1);
      my ($ss2,$mfe2) = RNA::fold($seq2);
      $D[$dim*($i-1)+($j-1)] =  $D[$dim*($j-1)+($i-1)] = RNA::bp_distance($ss1,$ss2);
      # printf "%s\n%s [ %6.2f ]\n", $seq1, $ss1, $mfe1;
      # printf "%s\n%s [ %6.2f ]\n", $seq2, $ss2, $mfe2;
      # print "bp_didtance $D[$dim*($i-1)+($j-1)]\n";
    }
    else {croak "Method $method not available ..please use SCI|dHn|dHx|dBp"}
  }

  # write matrix to file
  print STDERR "Writing distance matrix to file ...\n";
  if ($m eq "SCI"){ $Dfile = dump_matrix(\@D,$dim,1,1,"S",$ASO)}
  elsif ($m eq "dHn"){$Dfile = dump_matrix(\@D,$dim,1,1,"dHn",$ASO)}
  elsif ($m eq "dHx"){$Dfile = dump_matrix(\@D,$dim,1,1,"dHx",$ASO)}
  elsif ($m eq "dBp"){$Dfile = dump_matrix(\@D,$dim,1,1,"dBp",$ASO)}
  else { croak "Method $m not available ..please use SCI|dHn|dHx|dBp"}

  # check triangle inequality
  if ($check == 1){
    print STDERR "Checking triangle inequality ...\n";
    my $result = check_triangle($dim,\@D);
  }
  return $Dfile;
}


sub dump_matrix {
  my ($M,$d,$ad,$pd,$what,$aso) = @_;
  my ($i,$j,$info);
  my $ad_mx = file($outdir,"ld.mx"); # AnalyseDists lower diagoinal distance matrix
  my $pd_mx = file($outdir,"phylip.dst"); # Phylip distance matrix

  if ($what eq "S"){$info="> S (SCI distance)"}
  elsif($what eq "dHn"){$info="> H (Hamming distance with gap Ns)"}
  elsif($what eq "dHx"){$info="> H (Hamming distance with gaps removed)"}
  elsif($what eq "dBp"){$info="> B (Base pair distance)"}
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
      my $val = ${$aso->next_aln}{_order}->{$i};
      $val=~s/\//_/g;
      my $id = join '_', eval($i+1),$val;
      print $matrix $id." ";
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

A neighbour joining tree is reconstructed from pairwise distances of
sequences in the input alignment and subsets of the alignment are
derived by splitting at each edge of the NJ tree as well as performing
a split decomposition of the matrix of pairwise distances. These
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
'dHx', 'dBp', and 'SCI'. The first and second compute pairwise Hamming
distances of sequences, where 'dHn' replaces gaps with 'N', whereas
'dHx' removes all gap columns (is not yet implemented). 'dBp' folds
RNA sequences into thir MFE structures and computes pairwise base pair
distances. 'SCI' computes the distance as 1-log(SCI), based on a
truncated strucure conservation index of two sequences. The latter,
however, is not a metric and therefore often results in negative
branch lengths in Neighbor Joining trees. Use with caution. [default:
'dHn']

=item B<--out|-o>

Output base directory. Temporary data and results will be written to
this directory

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut
(END)

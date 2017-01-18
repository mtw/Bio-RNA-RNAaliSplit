#!/usr/bin/env perl
# Last changed Time-stamp: <2017-01-18 19:58:41 mtw>
# -*-CPerl-*-
#
# usage: alisplit.pl myfile.aln

use AlignSplit;
use WrapRNAz;
use SplitDecomposition;

use Data::Dumper;
use Carp;

my $alnfile = "./ALL.SL.short2.JEVG.1.2.mloc.aln";
my $format = "ClustalW";
# Display ID handling in Bio::AlignIO is broken for Stockholm format
# use ClustalW format instead !!!

my @nseqs=();
my ($i,$j,$dim,$Dfile);
my @pw_alns = ();
my @D = ();
my $method = "dHn"; # SCI | dHn | dHx

my $AlignSplitObject = AlignSplit->new(infile => $alnfile,
				       format => $format,
				       dump => 1);

#print Dumper($AlignSplitObject);
#die;

# extract all pairwise alignments
my $dim = $AlignSplitObject->next_aln->num_sequences;
#print "NS: $dim\n";
for ($i=1;$i<$dim;$i++){
  for($j=$i+1;$j<=$dim;$j++){
    my $sa = $AlignSplitObject->dump_subalignment("pairwise", undef, [$i,$j]);
    push @pw_alns, $sa->stringify;
  }
}

# initialize distance matrix
for($i=0;$i<$dim;$i++){
  for ($j=0;$j<$dim;$j++){
    $D[$dim*$i+$j] = 0.;
  }
}

# build distance matrix based on pairwise alignments
foreach my $ali (@pw_alns){
  #print "processing $ali ...\n";
  my $pw_aso = AlignSplit->new(infile => $ali,
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
    $Dfile = dump_matrix(\@D,$dim,1,"S");
  }
  elsif ($method eq "dHn") {
    # Hamming dist with gaps replaced by Ns
    $D[$dim*($i-1)+($j-1)] =  $D[$dim*($j-1)+($i-1)] = $pw_aso->hammingdistN;
    $Dfile = dump_matrix(\@D,$dim,1,"dHn");
  }
  elsif ($method eq "dHx") {
    # TODO compute $pw_aso->hammingdistX in AlignSplit.pm
    # Hamming dist with all gap columns removed
    $D[$dim*($i-1)+($j-1)] =  $D[$dim*($j-1)+($i-1)] = $pw_aso->hammingdistX;
    $Dfile = dump_matrix(\@D,$dim,1,"dHx");
  }
  else {
    croak "Method $method not available ..please use SCI|dHn|dHx";
  }
}

# compute Neighbor Joining tree and do split decomposition
my $sd = SplitDecomposition->new(infile => $Dfile,
				 odir => $AlignSplitObject->odir,
				 basename => $AlignSplitObject->infilebasename);
#print Dumper($sd);
print "Identified ".$sd->count." splits\n";
my $splitnr=1;
while (my $sets = $sd->pop()){
  my ($rnazo1,$rnazo2,$have_rnazo1,$have_rnazo2) = (0)x4;
  my $set1 = $$sets{S1};
  my $set2 = $$sets{S2};
  my $token = "split".$splitnr;
#  print "set1: @$set1\n";
#  print "set2: @$set2\n";
  $sa1 = $AlignSplitObject->dump_subalignment("splits", $token, $set1);
  $sa2 = $AlignSplitObject->dump_subalignment("splits", $token, $set2);
  if( scalar(@$set1) > 1){
    $rnazo1 = WrapRNAz->new(alnfile => $sa1);
    $have_rnazo1 = 1;
    #print Dumper($rnazo1);
    print join "\t",$rnazo1->P,$rnazo1->sci,$rnazo1->z,$sa1."\n";
  }
  if( scalar(@$set2) > 1){
    $rnazo2 = WrapRNAz->new(alnfile => $sa2);
    $have_rnazo2 = 1;
    print join "\t",$rnazo2->P,$rnazo2->sci,$rnazo2->z,$sa2."\n";
   # print "$rnazo2->P\t$rnazo2->sci\t$rnazo2->z\t$sa2\n";
    #print Dumper($rnazo2);
  }
  #print "----\n";
  $splitnr++;
}


###############
# subroutines #
###############


sub dump_matrix {
  my ($M, $d, $ad, $what) = @_;
  my ($i,$j);
  my $mxfile = "ld.mx";

  if ($what eq "S"){$info="> S (SCI distance)"}
  elsif($what eq "dHn"){$info="> H (Hamming distance with gap Ns)"}
  elsif($what eq "dHx"){$info="> H (Hamming distance with gaps removed)"}
  else{$info="> H (unknown)"}

  if (defined $ad){ # print lower triangle matrix input for AanalyseDists
    open my $matrix, ">", $mxfile or die $!;
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
  return $mxfile;
}

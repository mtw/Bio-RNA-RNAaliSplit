#!/usr/bin/env perl
# Last changed Time-stamp: <2017-01-13 21:15:16 mtw>
# -*-CPerl-*-
#
# usage: alisplit.pl myfile.aln

use AlignSplit;
use SplitDecomposition;
use Data::Dumper;
use Carp;

my $alnfile = "./t1_no11.aln";
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
    my $sa = $AlignSplitObject->dump_subalignment("pairwise", [$i,$j]);
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
  #  print "processing $ali ...\n";
  my $pw_aso = AlignSplit->new(infile => $ali,
			       format => "ClustalW");
#  print Dumper($pw_aso);
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
print Dumper($sd);


###############
# subroutines #
###############

sub build_NJtree {
  my $file = shift;
  print Dumper($file);
  my $AScmd = "cat $file | AnalyseDists -Xn";
  print "$AScmd\n";
}

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

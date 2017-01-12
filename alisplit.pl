#!/usr/bin/env perl
# Last changed Time-stamp: <2017-01-12 14:05:49 mtw>
# -*-CPerl-*-
#
# usage: alisplit.pl myfile.{aln,stk}

use AlignSplit;
use Data::Dumper;

my $alnfile = "./result.aln";
my $stkfile = "./ZIKV_SL1.stk";

my $format = "Stockholm";
my @nseqs=();
my ($i,$j,$dim);
my @pw_alns = ();
my @D = ();

my $AlignSplitObject = AlignSplit->new(infile => $stkfile,
				       format => $format,
				       dump => 1);

#print Dumper($AlignSplitObject);
#die;

# extract all pairwise alignments
my $dim = $AlignSplitObject->next_aln->num_sequences;
print "NS: $dim\n";
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

# build distance matrix based on -log( [normalized] pairwise SCI)
foreach my $ali (@pw_alns){
  my ($sci,$dsci);
  print "processing $ali ...\n";
  my $pw_aso = AlignSplit->new(infile => $ali,
			       format => "ClustalW");
#  my $bn = basename($pw_aso->infile->basename, ".aln");
  my ($i,$j) = sort split /_/, $pw_aso->infilebasename;

  if($pw_aso->sci > 1){$sci = 1}
  elsif ($pw_aso->sci == 0.){$sci = 0.000001}
  else { $sci = $pw_aso->sci; }
  $dsci = -1*log($sci)/log(10);

  $D[$dim*($i-1)+($j-1)] =  $D[$dim*($j-1)+($i-1)] = $dsci;
  print "$i -> $j : $sci\t$dsci\n";
}
dump_matrix(\@D,$dim,1);

sub dump_matrix {
  my ($M, $d, $ad) = @_;
  my ($i,$j);

  if (defined $ad){ # print lower triangle matrix input for AanalyseDists
    open my $matrix, ">", "ld.mx" or die $!;
    print $matrix "> S (SCI distance)\n";
    print $matrix "> X $d\n";
    for ($i=1;$i<$d;$i++){
      for($j=0;$j<$i;$j++){
	printf $matrix "%6.4f ", @$M[$d*$i+$j];
      }
      print $matrix "\n";
    }
    close $matrix;
  }
  else{ # print the entire matrix
    for ($i=0;$i<$d;$i++){
      for($j=0;$j<$d;$j++){
	printf "%6.4f ", @$M[$d*$i+$j];
      }
      print "\n";
    }
  }
  
}

die;










my $aln = $AlignSplitObject->alignment->next_aln();
# Extract sequences and check values for the alignment column $pos
foreach $seq ($aln->each_seq) {
  $_ = $seq->seq;
  s/-/N/g;
  push @nseqs, $_;
}

my $allseq = join "\n", @nseqs;
my $AnalyseSequences_cmd = "echo \'$allseq\' | AnalyseSeqs -Xm";
#print $AnalyseSequences_cmd."\n";

open (AS, $AnalyseSequences_cmd."|");
while(<AS>){
  chomp;
  next if($.==1);
  if($.==2){
    print "\@\@$_\@\@\n";
    die unless ($_ =~ /^>\sX\s+(\d+)/);
    $dim = $1;
    print "dim is $dim\n";
    next
  }
  print $.;
}
close(AS);
#### TODO read in distance matrix produced by another program, such as AnalyseSequences

#my $parser = Bio::Matrix::IO->new(-format => 'phylip',
#                                    -file   => 'filename.dist');
#my $mat  = $parser->next_matrix;
#my $tree = $dfactory->make_tree($mat);
#$treeout->write_tree($tree);

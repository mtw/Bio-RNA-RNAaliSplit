# -*-CPerl-*-
# Last changed Time-stamp: <2020-02-16 05:52:03 mtw>
#
#  Derive features of an alignment, in particular scores to compare
#  different alignments of the same sequences

package Bio::RNA::RNAaliSplit::AliFeature;

use Moose;
use namespace::autoclean;
use version; our $VERSION = qv('0.11');
use warnings;
use Data::Dumper;
use Storable 'dclone';
use Carp;
use RNA;

extends 'Bio::RNA::RNAaliSplit::AliHandler';

has 'sop' => ( # sum of pairs score
	      is => 'rw',
	      isa => 'Int',
	      predicate => 'has_sSOP',
	      init_arg => undef,
	     );

has '_csp' => ( # column sequence positions
	      is => 'ro', # read-only
	      isa => 'ArrayRef',
	      predicate => 'has_sCSP',
	      init_arg => undef,
	      writer => '_cspwriter', # private writer
	      );

has 'csp_hash' => (
		   is => 'ro',
		   isa => 'HashRef',
		   predicate => 'hash_csp_hash',
		   init_arg => undef,
		   writer => '_csp_hash_writer', # private writer 4 ro attribute
		  );

has 'hammingdistN' => (
		       is => 'rw',
		       isa => 'Num',
		       default => '-1',
		       predicate => 'has_hammingN',
		       init_arg => undef,
		      );

has 'hammingdistX' => (
		       is => 'rw',
		       isa => 'Num',
		       default => '-1',
		       predicate => 'has_hammingX',
		       init_arg => undef,
		      );

has 'sci' => (
              is => 'rw',
              isa => 'Num',
              init_arg => undef,
              documentation => q(Structure conservation index),
             );

with 'FileDirUtil';

sub BUILD {
  my $self = shift;
   my $this_function = (caller(0))[3];
  confess "ERROR [$this_function] \$self->ifile not available"
    unless ($self->has_ifile);
  $self->alignment({-file => $self->ifile,
		    -format => $self->format,
		    -displayname_flat => 1} ); # discard position in sequence IDs
  $self->next_aln($self->alignment->next_aln);
  $self->next_aln->set_displayname_safe();
  $self->_get_alen();
  $self->_get_nrseq();
  $self->set_ifilebn;
  if ($self->nrseq == 2){ $self->_hamming() }

  # compute sum of pairs score
  #$self->compute_sop();
  # compute sequence position for each column
  $self->_get_column_sequence_positions();
  # compute CSP hash
  $self->_csp_hash();
	$self->_compute_sci();
}

# Compute Sum of Pairs scoring for an alignment
sub compute_sop {
  my $self = shift;
  # print "### in ccompute_sop###\n";
  my @alignment=();
  my $sp=0;
  for (my $i=1;$i<=$self->nrseq;$i++){
    push @alignment, $self->next_aln->get_seq_by_pos($i)->seq;
  }
  # print Dumper (\@alignment);

  for (my $c=0;$c<eval($self->alen);$c++){ # loop over alignment columns
    my $colscore=0;
    # print "$c\n";

    for (my $i=0;$i<eval($self->nrseq);$i++){
      for (my $j=$i+1;$j<eval($self->nrseq);$j++){
 	my $ci = substr($alignment[$i],$c,1);
	my $cj = substr($alignment[$j],$c,1);
	if (1) { # if ($self->distancemeasure eq "e"){ # edit distance
 	  if ($ci ne $cj) {$colscore += 1}
 	}
 	# print " > $ci$i $cj$j $colscore <\n";
      } # end for j
    } # end for i
    # print " > ---------- <\n";
    $sp += $colscore;
  } # end for c
  $self->sop($sp);
}


sub _get_column_sequence_positions {
  my $self = shift;
  my @loclist = ();
  foreach my $i( 1..$self->nrseq ) {
    my @ll=();
    my $seq = $self->next_aln->get_seq_by_pos($i);
    # print $seq->seq."\n";
    $ll[0] =  $self->alen;
    for (my $j=1;$j<=$self->alen;$j++){
      my $pos=0; #default
      my $loc = $seq->location_from_column($j);
      if (defined ($loc)){
	if($loc->location_type() eq 'EXACT'){
	  $pos = $loc->to_FTstring();
	}
	elsif ($loc->location_type() eq 'IN-BETWEEN'){
	  $pos = 0;
	}
	else { croak "ERROR: this should not happen\n".Dumper($loc); }
      }
      else {
	  $pos = 0; # TODO check me
      }
      # print Dumper($loc);
      $ll[$j]=$pos;
    } # end for
    push @loclist, \@ll;
  } # end foreach
  $self->_cspwriter(\@loclist);
}

sub _csp_hash {
  my $self = shift;
  my %csp = ();
  for (my $j=1;$j<=$self->alen;$j++) { # loop over columns
    my $pstring;
    for (my $i=0;$i<$self->nrseq;$i++){ # loop over sequences
      $pstring .= eval(${$self->_csp}[$i]->[$j]).":";
    }
    #print ">> $pstring <<\n";
    unless (defined $csp{$pstring}){ $csp{$pstring}=1;}
    else {$csp{$pstring}+=1;}
  }
  $self->_csp_hash_writer(\%csp);
}

sub _hamming {
  my $self = shift;
  my $this_function = (caller(0))[3];
  my $hamming = -1;
  croak "ERROR [$this_function] cannot compute Hamming distance for $self->next_aln->num_sequences sequences"
    if ($self->next_aln->num_sequences != 2);

  my $aln =  $self->next_aln->select_noncont((1,2));

  # compute Hamming distance of the aligned sequences, replacing gaps with Ns
  my $alnN = dclone($aln);
  croak("ERROR [$this_function] cannot replace gaps with Ns")
    unless ($alnN->map_chars('-','N') == 1);
  my $seq1 = $alnN->get_seq_by_pos(1)->seq;
  my $seq2 = $alnN->get_seq_by_pos(2)->seq;
  croak "ERROR [$this_function] sequence length differs"
    unless(length($seq1)==length($seq2));
  my $hammingN = ($seq1 ^ $seq2) =~ tr/\001-\255//;
  $self->hammingdistN($hammingN);
}

sub _compute_sci {
  my $self = shift;
	my $this_function = (caller(0))[3];
  my ($aln,$sci, $amfe, $avmfe, $afc, $cs);
	my @aln=();
  $aln = $self->next_aln();
  my $nrseq = 0;
  my $sum = 0;
  # Extract sequences and check values for the alignment column $pos
  foreach my $seqs ($aln->each_seq) {
    my $seq = $seqs->seq;
    my $id = $seqs->id;
    push @aln,$seq;
    my $fc= new RNA::fold_compound($seq); # compute invididual MFEs
    my ($ss,$mfe) = $fc->mfe();
    #print ">$id\n$seq\n$ss ($mfe)\n---\n";
    $nrseq++;
    $sum += $mfe;
  }
  #print "$nrseq sequences\n";
  $avmfe = $sum/$nrseq;
  $afc =  new RNA::fold_compound(\@aln);
  ($cs,$amfe) = $afc->mfe(); #alifold consensus mfe
  $self->sci($amfe/$avmfe);
  #print "selfcomputed SCI: ".$self->sci."\n";# ($amfe / $avmfe)\n";
  return ($cs,$sci);
}


no Moose;

1;

__END__

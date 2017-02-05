# -*-CPerl-*-
# Last changed Time-stamp: <2017-02-05 23:23:26 mtw>

# AlignSplit.pm: Handler for horizontally splitting alignments
#
# requires RNAalifold executable available to the perl interpreter

package AlignSplit;

use version; our $VERSION = qv('0.02');
use Carp;
use Data::Dumper;
use Moose;
use Moose::Util::TypeConstraints;
use Path::Class::File;
use Path::Class::Dir;
use Path::Class;
use File::Basename;
use IPC::Cmd qw(can_run run);
use Bio::AlignIO;
use Storable 'dclone';

subtype 'MyAln' => as class_type('Bio::AlignIO');

coerce 'MyAln'
    => from 'HashRef'
    => via { Bio::AlignIO->new( %{ $_ } ) };

has 'format' => ( 
		 is => 'ro',
		 isa => 'Str',
		 predicate => 'has_format',
		 default => 'ClustalW',
		 required => 1,
		);

has 'infilebasename' => (
			 is => 'rw',
			 isa => 'Str',
			 predicate => 'has_basename',
			 init_arg => undef,
			);

has 'alignment' => (
		    is => 'rw',
		    isa => 'MyAln',
		    predicate => 'has_aln',
		    coerce => 1,
		   );

has 'next_aln' => (
		   is => 'rw',
		   isa => 'Bio::SimpleAlign',
		   predicate => 'has_next_aln',
		   init_arg => undef,
		  );

has 'dump' => (
	       is => 'rw',
	       isa => 'Num',
	       predicate => 'has_dump_flag',
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

with 'FileDir';

sub BUILD {
    my $self = shift;
    my $this_function = (caller(0))[3];
    confess "ERROR [$this_function] \$self->ifile not available"
      unless ($self->has_ifile);
    $self->alignment({-file => $self->ifile,
		      -format => $self->format,
		      -displayname_flat => 1} ); # discard position in sequence IDs
    $self->next_aln($self->alignment->next_aln);
    $self->odir( [$self->ifile->dir,$self->odirn] );
    mkdir($self->odir);
    $self->infilebasename(fileparse($self->ifile->basename, qr/\.[^.]*/));

    if ($self->has_dump_flag){
      # dump ifile as aln in ClustalW format to odir/input
      my $iodir = $self->odir->subdir('input');
      mkdir($iodir);
      my $ialnfile = file($iodir,$self->infilebasename.".aln");
      my $alnio = Bio::AlignIO->new(-file   => ">$ialnfile",
				    -format => "ClustalW",
				    -flush  => 0,
				    -displayname_flat => 1 );
      my $aln2 = $self->next_aln->select_noncont((1..$self->next_aln->num_sequences));
      $alnio->write_aln($aln2);
      # end dump input aln file
    }

    if ($self->next_aln->num_sequences == 2){ $self->_hamming() }
  }

sub dump_subalignment {
  my ($self,$alipathsegment,$token,$what) = @_;
  my $this_function = (caller(0))[3];
  my ($aln,$aln2,$name);

  croak "ERROR [$this_function] argument 'token' not provided"
    unless (defined($token));

  # create output path
  my $ids = join "_", @$what;
  unless (defined($alipathsegment)){$alipathsegment = "tmp"}
  my $oodir = $self->odir->subdir($alipathsegment);
  mkdir($oodir);

  # create info file 
  my $oinfofile = file($oodir,$token.".info");
  open my $oinfo, ">", $oinfofile or die $!;
  foreach my $entry (@$what){
    my $key = $entry-1;
    my $val = ${$self->next_aln}{_order}->{$key};
    print $oinfo join "\t", ($entry, $val, "\n");
  }
  close($oinfo);

  # create subalignment in Clustal and Stockholm format
  my $oalifile_clustal = file($oodir,$token.".aln");
  my $oalifile_stockholm = file($oodir,$token.".stk");
  my $oali_clustal = Bio::AlignIO->new(-file   => ">$oalifile_clustal",
				       -format => "ClustalW",
				       -flush  => 0,
				       -displayname_flat => 1 );
  my $oali_stockholm = Bio::AlignIO->new(-file   => ">$oalifile_stockholm",
					 -format => "Stockholm",
					 -flush  => 0,
					 -displayname_flat => 1 );
  $aln = $self->next_aln->select_noncont(@$what);
  $oali_clustal->write_aln( $aln );
  $oali_stockholm->write_aln( $aln );

  # create subalignment fasta file
  my $ofafile = file($oodir,$token.".fa");
  my $ofa = Bio::AlignIO->new(-file   => ">$ofafile",
			      -format => "fasta",
			      -flush  => 0,
			      -displayname_flat => 1 );
  $aln2 = $aln->remove_gaps;
  $ofa->write_aln( $aln2 );

  # extract sequences from alignment and dump to .seq file
  # NOTE that these sequences do contain gap symbols, intentionally (!)
  # these can then be replaced by Ns to compute eg hamming distance
  my $oseqfile = file($oodir,$token.".seq");
  open my $seqfile, ">", $oseqfile or die $!;
  foreach my $seq ($aln->each_seq) {
    print $seqfile $seq->seq,"\n";
  }
  close($seqfile);

  return ( $oalifile_clustal,$oalifile_stockholm );
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

#  print $self->infilebasename,":\n";
#  print ">>s1: $seq1\n";
#  print ">>s2: $seq2\n";
#  print "** dhN = ".$self->hammingdistN."\n";
#  print "+++\n";
}

no Moose;

1;

# -*-CPerl-*-
# Last changed Time-stamp: <2019-04-10 09:56:57 mtw>
#
#  Derive features of an alignment, in particular scores to compare
#  different alignments of the same sequences

package Bio::RNA::RNAaliSplit::AliFeature;

use Moose;
use namespace::autoclean;
use version; our $VERSION = qv('0.11');
use diagnostics;
use Data::Dumper;

extends 'Bio::RNA::RNAaliSplit::AliHandler';

has 'sop' => ( # sum of pairs score
	      is => 'rw',
	      isa => 'Int',
	      predicate => 'has_sSOP',
	      init_arg => undef,
	     );

has 'cols' => ( # column score
	       is => 'rw',
	       isa => 'Arrayref',
	       predicate => 'has_sCOL',
	       init_arg => undef,
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
  # ev past odir stuff here
  $self->set_ifilebn;

  # compute sum of pairs score
  #$self->sop();
  # compute column score
  #$self->cols();
}

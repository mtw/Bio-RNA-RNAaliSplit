# -*-CPerl-*-
# Last changed Time-stamp: <2016-12-31 22:42:05 mtw>
package AlignSplit;

use version; our $VERSION = qv('0.01');
use Carp;
use Data::Dumper;
use Moose;
use Moose::Util::TypeConstraints;
use Path::Class::File;
use Bio::AlignIO;

subtype 'MyAln' => as class_type('Bio::AlignIO');

coerce 'MyAln'
    => from 'HashRef'
    => via { Bio::AlignIO->new( %{ $_ } ) };

subtype 'MyFile' => as class_type('Path::Class::File');

coerce 'MyFile'
    => from 'Str'
    => via { Path::Class::File->new($_) };

has 'file' => (
    is => 'ro',
    isa => 'MyFile',
    predicate => 'has_file',
 #   required => 1,
    coerce => 1,
    );

has 'format' => ( 
    is => 'ro',
    isa => 'Str',
    predicate => 'has_format',
    default => 'ClustalW',
    required => 0,
    );

has 'alignment' => (
    is => 'rw',
    isa => 'MyAln',
    predicate => 'has_aln',
    coerce => 1,
    );

has 'basename' => (
    is => 'rw',
    isa => 'Str',
    predicate => 'has_basename',
    );

has 'path' => (
    is => 'rw',
    isa => 'Str',
    predicate => 'has_path',
    );

has 'ext' => (
    is => 'rw',
    isa => 'Str',
    predicate => 'has_ext',
    );


sub BUILD {
    my $self = shift;
    my $this_function = (caller(0))[3];
  #  my ($x,$y,$z);
    carp "INFO [$this_function]";
    
    confess "ERROR [$this_function] \$self->file not available"
	unless ($self->has_file);
 
    $self->alignment({
	-file => $self->file,
	-format => $self->format});
    
  #  ($x,$y,$z) = fileparse($self->bamfile,qr /\..*/);
  #  $self->basename($x);
  #  $self->path($y);
  #  $self->ext($z);
}


no Moose;

1;

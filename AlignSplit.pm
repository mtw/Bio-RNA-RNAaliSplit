# -*-CPerl-*-
# Last changed Time-stamp: <2017-01-10 10:09:33 mtw>
package AlignSplit;

use version; our $VERSION = qv('0.01');
use Carp;
use Data::Dumper;
use Moose;
use Moose::Util::TypeConstraints;
use Path::Class::File;
use IPC::Cmd qw(can_run run);
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

has 'sci' => (
	      is => 'rw',
	      isa => 'Num',
	      predicate => 'has_sci',
	      builder => '_compute_sci',
	      lazy => 1, # required to make sure $self->file is available
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

sub _compute_sci {
  my $self = shift;
  my $this_function = (caller(0))[3];
  my $sci = 0.0;
  my $f;
  my $alifold = can_run('RNAalifold') or
    croak " ERROR [$this_function] RNAalifold not found";
 # carp "[$this_function]: $alifold";

  $f = $self->file->stringify;
  my $cmd = "$alifold --sci $f";

  my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmd, verbose => 0 );
  if( !$success ) {
    print STDERR "ERROR [$this_function] Call to $alifold unsuccessful\n";
    print STDERR "ERROR: this is what the command printed:\n";
    print join "", @$full_buf;
    croak $!;
  }
  my @aliout = split '\n', $$stdout_buf[0];
  print Dumper(@aliout);
  print ">>$aliout<<\n";
  print Dumper($full_buf);
  return $sci
}


no Moose;

1;

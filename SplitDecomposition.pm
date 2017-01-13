# -*-CPerl-*-
# Last changed Time-stamp: <2017-01-13 21:14:46 mtw>

# SplitDecomposition.pm: Wrapper for computing split decomposition
#
# requires AnalyzeDists executable from the ViennaRNA poackage
# available to the perl interpreter

package SplitDecomposition;

use version; our $VERSION = qv('0.01');
use Carp;
use Data::Dumper;
use Moose;
use Moose::Util::TypeConstraints;
use Path::Class::File;
use Path::Class::Dir;
use Path::Class;
#use File::Basename;
use IPC::Cmd qw(can_run run);

my ($analysedists,$oodir);

subtype 'MyFile2' => as class_type('Path::Class::File');

coerce 'MyFile2'
  => from 'Str'
  => via { Path::Class::File->new($_) };

subtype 'MyDir2' => as class_type('Path::Class::Dir');

coerce 'MyDir2'
  => from 'ArrayRef'
  => via { Path::Class::Dir->new( @{ $_ } ) };

has 'infile' => (
		 is => 'ro',
		 isa => 'MyFile2',
		 predicate => 'has_infile',
		 coerce => 1,
		 required => 1,
		 documentation => q(lower triangular matrix of distances),
		);

has 'basename' => (
		   is => 'rw',
		   isa => 'Str',
		   predicate => 'has_basename',
		   );

has 'odirname' => (
		   is => 'ro',
		   default => 'as',
		   predicate => 'has_odirname',
		  );

has 'odir' => (
	       is => 'rw',
	       isa => 'MyDir2',
	       predicate => 'has_odir',
	       coerce => 1,
	      );

sub BUILD {
  my $self = shift;
  my $this_function = (caller(0))[3];
  confess "ERROR [$this_function] \$self->infile not available"
    unless ($self->has_infile);
   $analysedists = can_run('AnalyseDists') or
    croak " ERROR [$this_function] AnalyseDists not found";
  unless($self->has_odir){
    unless($self->has_odirname){
      self->odirname("as");
    }
    $self->odir( [$self->infile->dir,$self->odirname] );
    mkdir($self->odir);
  }
  $oodir = $self->odir->subdir("analysedists");
  mkdir($oodir);

  # do computation
  $self->_NeighborJoining();
}

sub _NeighborJoining {
  # TODO  warn if negative branch lengths occur`
  my $self = shift;
  my $this_function = (caller(0))[3];
  my ($nj_outfilename, $nj_treefilename);

  if ($self->has_basename){
    $nj_outfilename = $self->basename.".nj.out";
    $nj_treefilename = $self->basename.".nj.ps";
  }
  else{
    $nj_outfilename = "nj.out";
    $nj_treefilename = "nj.ps";
  }
  my $nj_out = file($oodir,$nj_outfilename);
  my $nj_tree = file($oodir,$nj_treefilename);
  open my $fh, ">", $nj_out;

  my $ad_cmd = $analysedists." -Xn < ".$self->infile;
  my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $ad_cmd, verbose => 0 );
  if( !$success ) {
    print STDERR "ERROR [$this_function] Call to $analysedists unsuccessful\n";
    print STDERR "ERROR: this is what the command printed:\n";
    print join "", @$full_buf;
    croak $!;
  }
  my @adout = split /\n/, $$stdout_buf[0];
  foreach my $line( @adout){
    print $fh $line,"\n";
  }
  close($fh);
  rename "nj.ps", $nj_tree;
}

1;


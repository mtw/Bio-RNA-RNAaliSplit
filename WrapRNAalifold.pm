# -*-CPerl-*-
# Last changed Time-stamp: <2017-01-23 19:22:28 mtw>

# WrapRNAz.pm: A versatile object-oriented wrapper for RNAalifold
#
# Requires RNAalifold executable from the ViennaRNA package available
# to the Perl interpreter.

package WrapRNAalifold;

use version; our $VERSION = qv('0.01');
use Carp;
use Data::Dumper;
use Moose;
use Moose::Util::TypeConstraints;
use Path::Class::File;
use Path::Class::Dir;
use Path::Class;
use File::Basename;
use IPC::Cmd qw(can_run run);

my ($rnaalifold,$oodir);

subtype 'MyFile4' => as class_type('Path::Class::File');

coerce 'MyFile4'
  => from 'Str'
  => via { Path::Class::File->new($_) };

subtype 'MyDir4' => as class_type('Path::Class::Dir');

coerce 'MyDir4'
  => from 'ArrayRef'
  => via { Path::Class::Dir->new( @{ $_ } ) };

has 'alnfile' => (
		 is => 'ro',
		 isa => 'MyFile4',
		 predicate => 'has_alnfile',
		 coerce => 1,
		 required => 1,
		 documentation => q(Alignment file in ClustalW format),
		 );

has 'alnfilebasename' => (
			  is => 'rw',
			  isa => 'Str',
			  predicate => 'has_alnfilebasename',
			  init_arg => undef, # make this unsettable via constructor
			 );

has 'bn' => (
	     is => 'rw',
	     isa => 'Str',
	     predicate => 'has_basename',
	     documentation => q(Set this to override output basename),
	    );

has 'odirname' => (
		   is => 'ro',
		   default => 'results',
		   predicate => 'has_odirname',
		  );

has 'odir' => (
	       is => 'rw',
	       isa => 'MyDir4',
	       predicate => 'has_odir',
	       coerce => 1,
	      );

has 'P' => (
	    is => 'rw',
	    isa => 'Num',
	    init_arg => undef,
	    documentation => q(SVM RNA-class probability),
	   );

has 'z' => (
	    is => 'rw',
	    isa => 'Num',
	    init_arg => undef,
	    documentation => q(Mean z-score),
	   );

has 'sci' => (
	      is => 'rw',
	      isa => 'Num',
	      init_arg => undef,
	      documentation => q(Structure conservation index),
	     );

has 'consensus_struc' => (
			  is => 'rw',
			  isa => 'Str',
			  predicate => 'has_consensus_struc',
			  init_arg => undef,
			 );

has 'consensus_mfe' => (
			is => 'rw',
			isa => 'Num',
			predicate => 'has_consensus_mfe',
			init_arg => undef,
		       );

has 'consensus_energy' => (
			   is => 'rw',
			   isa => 'Num',
			   predicate => 'has_consensus_energy',
			   init_arg => undef,
			  );

has 'consensus_covar_terms' => (
				is => 'rw',
				isa => 'Num',
				predicate => 'has_consensus_covar_terms',
				init_arg => undef,
			       );


sub BUILD {
  my $self = shift;
  my $this_function = (caller(0))[3];
  confess "ERROR [$this_function] \$se+lf->alnfile not available"
    unless ($self->has_alnfile);
   $rnaalifold = can_run('RNAalifold') or
     croak "ERROR [$this_function] RNAalifold not found";
  unless($self->has_odir){
    $self->odir( [$self->alnfile->dir,$self->odirname] );
    mkdir($self->odir);
  }
  $oodir = $self->odir->subdir("alifold");
  mkdir($oodir);
  $self->alnfilebasename(fileparse($self->alnfile->basename, qr/\.[^.]*/));

  $self->run_rnaalifold();
}

sub run_rnaalifold {
  my $self = shift;
  my $this_function = (caller(0))[3];
  my ($out_fn,$out,$alnps_fn,$alnps,$alirnaps_fn,$alirnaps,$alidotps_fn,$alidotps);
  if ($self->has_alnfilebasename){
    $out_fn = $self->alnfilebasename.".alifold.out";
    $alnps_fn = $self->alnfilebasename.".aln.ps";
    $alirnaps_fn = $self->alnfilebasename.".alirna.ps";
    $alidotps_fn = $self->alnfilebasename.".alidot.ps";
  }
  elsif ($self->has_basename){
    $out_fn = $self->bn.".alifold.out";
    $alnps_fn = $self->bn.".aln.ps";
    $alirnaps_fn = $self->bn.".alirna.ps";
    $alidotps_fn = $self->bn.".alidot.ps";
  }
  else{
    $out_fn = "alifold.out";
    $alnps_fn = "aln.ps";
    $alirnaps_fn = "alirna.ps";
    $alidotps_fn = "alidot.ps";
  }
  $out = file($oodir,$out_fn); # RNAalifold stdout
  $alnps = file($oodir,$alnps_fn); # RNAalifold aln.ps
  $alirnaps = file($oodir,$alirnaps_fn); # RNAalifold alirna.ps
  $alidotps = file($oodir,$alidotps_fn); # RNAalifold alidot.ps

  open my $fh, ">", $out;
  my $cmd = $rnaalifold." --aln --color -r --cfactor 0.6 --nfactor 0.5 -p --sci ".$self->alnfile;
  my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmd, verbose => 0 );
  if( !$success ) {
    print STDERR "ERROR [$this_function] Call to $rnaalifold unsuccessful\n";
    print STDERR "ERROR: this is what the command printed:\n";
    print join "", @$full_buf;
    croak $!;
  }
  my $stdout_buffer = join "", @$stdout_buf;
  my @rnaalifoldout = split /\n/, $stdout_buffer;
  foreach my $line( @rnaalifoldout){
    print $fh $line,"\n";
  }
  close($fh);

  $self->_parse_rnaalifold($stdout_buffer);
  rename "aln.ps", $alnps;
  rename "alirna.ps", $alirnaps;
  rename "alidot.ps", $alidotps;
  unlink "alifold.out";
}

# parse RNAalifold output
sub _parse_rnaalifold {
  my ($self,$out) = @_;
  my $this_function = (caller(0))[3];
  my @buffer=split(/^/, $out);

  foreach my $i (0..$#buffer){
    next unless ($i == 1); # just parse consensus structure
    unless ($buffer[$i] =~ m/([\(\)\.]+)\s+\(\s?(-?\d+\.\d+)\s+=\s+(-?\d+\.\d+)\s+\+\s+(-?\d+\.\d+)\)\s+\[sci\s+=\s+(\d+\.\d+)\]/){
      carp "ERROR [$this_function]  cannot parse RNAalifold output:";
      croak $self->alnfile.":\n$buffer[$i]";
    }
    $self->consensus_struc($1);
    $self->consensus_mfe($2);
    $self->consensus_energy($3);
    $self->consensus_covar_terms($4);
    $self->sci($5);
    last;
  }
}

1;

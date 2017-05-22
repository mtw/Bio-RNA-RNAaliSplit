# -*-CPerl-*-
# Last changed Time-stamp: <2017-05-15 12:24:42 mtw>

# Bio::RNA::RNAaliSplit::WrapRNAalifold.pm: A versatile object-oriented
# wrapper for RNAalifold
#
# Requires RNAalifold executable from the ViennaRNA package available
# to the Perl interpreter.

package Bio::RNA::RNAaliSplit::WrapRNAalifold;

use version; our $VERSION = qv('0.05_02');
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

has 'ribosum' => (
		  is => 'rw',
		  isa => 'Num',
		  predicate => 'has_ribosum',
		  documentation => q(Ribosum scoring),
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

has 'RNAalifold_version' => (
			     is => 'rw',
			     isa => 'String',
			     init_arg => undef,
			    );

with 'Bio::RNA::RNAaliSplit::FileDir';

sub BUILD {
  my $self = shift;
  my $this_function = (caller(0))[3];
  confess "ERROR [$this_function] \$self->ifile not available"
    unless ($self->has_ifile);
  $rnaalifold = can_run('RNAalifold') or
    croak "ERROR [$this_function] RNAalifold not found";
  unless($self->has_odir){
    unless($self->has_odirn){self->odirn("as")}
    $self->odir( [$self->ifile->dir,$self->odirn] );
    mkdir($self->odir);
  }
  $oodir = $self->odir->subdir("alifold");
  mkdir($oodir);
  $self->alnfilebasename(fileparse($self->ifile->basename, qr/\.[^.]*/));

  $self->run_rnaalifold();
}

sub run_rnaalifold {
  my $self = shift;
  my $this_function = (caller(0))[3];
  my ($out_fn,$out,$alnps_fn,$alnps,$alirnaps_fn,$stk_fn);
  my ($alirnaps,$alidotps_fn,$alidotps,$alifoldstk);
  my $tag = "";
  if ($self->has_ribosum){$tag = ".risobum"}
  if ($self->has_alnfilebasename){
    $out_fn = $self->alnfilebasename.$tag."."."alifold.out";
    $alnps_fn = $self->alnfilebasename.$tag."."."aln.ps";
    $alirnaps_fn = $self->alnfilebasename.$tag."."."alirna.ps";
    $alidotps_fn = $self->alnfilebasename.$tag."."."alidot.ps";
    $stk_fn = $self->alnfilebasename.$tag."."."alifold.stk";
  }
  elsif ($self->has_basename){
    $out_fn = $self->bn.$tag."."."alifold.out";
    $alnps_fn = $self->bn.$tag."."."aln.ps";
    $alirnaps_fn = $self->bn.$tag."."."alirna.ps";
    $alidotps_fn = $self->bn.$tag."."."alidot.ps";
    $stk_fn = $self->bn.$tag."."."alifold.stk";
  }
  else{
    $out_fn = $tag."alifold.out";
    $alnps_fn = $tag."aln.ps";
    $alirnaps_fn = $tag."alirna.ps";
    $alidotps_fn = $tag."alidot.ps";
  }
  $out = file($oodir,$out_fn); # RNAalifold stdout
  $alnps = file($oodir,$alnps_fn); # RNAalifold aln.ps
  $alirnaps = file($oodir,$alirnaps_fn); # RNAalifold alirna.ps
  $alidotps = file($oodir,$alidotps_fn); # RNAalifold alidot.ps
  $alifoldstk = file($oodir,$stk_fn); # RNAalifold generated Stockholm file with new CS 

  open my $fh, ">", $out;
  my $alifold_options = " --aln --color --cfactor 0.6 --nfactor 0.5 -p --sci --aln-stk ";
  if ($self->has_ribosum){$alifold_options.=" -r "}
  my $cmd = $rnaalifold.$alifold_options.$self->ifile;
  my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmd, verbose => 0 );
  if( !$success ) {
    print STDERR "ERROR [$this_function] Call to $rnaalifold unsuccessful\n";
    print STDERR "ERROR: $cmd\n";
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
  rename "RNAalifold_results.stk", $alifoldstk;
  unlink "alifold.out";
}

# parse RNAalifold output
sub _parse_rnaalifold {
  my ($self,$out) = @_;
  my $this_function = (caller(0))[3];
  my @buffer=split(/^/, $out);

  foreach my $i (0..$#buffer){
    next unless ($i == 1); # just parse consensus structure
    unless ($buffer[$i] =~ m/([\(\)\.]+)\s+\(\s*(-?\d+\.\d+)\s+=\s+(-?\d+\.\d+)\s+\+\s+(-?\d+\.\d+)\)\s+\[sci\s+=\s+(-?\d+\.\d+)\]/){
      carp "ERROR [$this_function]  cannot parse RNAalifold output:";
      croak $self->ifile.":\n$buffer[$i]";
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

# -*-CPerl-*-
# Last changed Time-stamp: <2017-03-11 17:20:44 michl>
# place of birth: somewhere over Newfoundland

# Bio::RNA::RNAaliSplit::WrapRscape.pm: A versatile object-oriented
# wrapper for R-scape
#
# Requires R-scape executable available to the Perl interpreter.

package Bio::RNA::RNAaliSplit::WrapRscape;

use version; our $VERSION = qv('0.04');
use Carp;
use Data::Dumper;
use Moose;
use Moose::Util::TypeConstraints;
use Path::Class::File;
use Path::Class::Dir;
use Path::Class;
use File::Basename;
use IPC::Cmd qw(can_run run);

my ($rscape,$oodir);

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

has 'statistic' => (
		is => 'rw',
		isa => 'String',
		predicate => 'has_statistic',
		documentation => q(Covariation statistic),
	     );

with 'Bio::RNA::RNAaliSplit::FileDir';

sub BUILD {
  my $self = shift;
  my $this_function = (caller(0))[3];
  confess "ERROR [$this_function] \$self->ifile not available"
    unless ($self->has_ifile);
  $rscape = can_run('R-scape') or
    croak "ERROR [$this_function] R-scape not found";
  unless($self->has_odir){
    unless($self->has_odirn){self->odirname("as")}
    $self->odir( [$self->ifile->dir,$self->odirn] );
    mkdir($self->odir);
  }
  $oodir = $self->odir->subdir("rscape");
  mkdir($oodir);
  $self->alnfilebasename(fileparse($self->ifile->basename, qr/\.[^.]*/));

  $self->run_rscape();
}

sub run_rscape {
  my $self = shift;
  my $this_function = (caller(0))[3];
  my ($bname,$out_fn,$sout_fn,$out,$sout,$sum_fn,$sum);
  my ($R2Rsto_fn,$R2Rsto,$R2Rstopdf_fn,$R2Rstopdf,$R2Rstosvg_fn,$R2Rstosvg);
  my $tag = "";
  if ($self->has_statistic){$tag = ".".$self->statistic};
  print ">> self->statistic is ".$self->statistic;die;

  if ($self->has_alnfilebasename){
    $bname = $self->alnfilebasename;
    $out_fn = $bname.$tag."."."rscape.out";
    $sout_fn = $bname.$tag."."."rscape.sorted.out";
    $sum_fn = $bname.$tag."."."rscape.sum";
    $R2Rsto_fn = $bname.$tag."."."R2R.sto";
    $R2Rstopdf_fn = $bname.$tag."."."R2R.sto.pdf";
    $R2Rstosvg_fn = $bname.$tag."."."R2R.sto.svg";
  }
  elsif ($self->has_basename){
    $bname = $self->bn;
    $out_fn = $bname.$tag."."."rscape.out";
    $sout_fn = $bname.$tag."."."rscape.sorted.out";
    $sum_fn = $bname.$tag."."."rscape.sum";
    $R2Rsto_fn = $bname.$tag."."."R2R.sto";
    $R2Rstopdf_fn = $bname.$tag."."."R2R.sto.pdf";
    $R2Rstosvg_fn = $bname.$tag."."."R2R.sto.svg";
  }
  else{
    $out_fn = $tag."rscape.out";
    $sout_fn = $tag."rscape.sorted.out";
    $sum_fn = $tag."rscape.sum";
    $R2Rsto_fn = $tag."R2R.sto";
    $R2Rstopdf_fn = $tag."R2R.sto.pdf";
    $R2Rstosvg_fn = $tag."R2R.sto.svg";
  }
  $out = file($oodir,$out_fn); # R-scape stdout
  $sout = file($oodir,$sout_fn); # R-scape sorted stdout
  $sum = file($oodir,$sum_fn); # R-scape .sum output
  $R2Rsto = file($oodir,$R2Rsto_fn); # R-scape R2R Stockholm file
  $R2Rstopdf = file($oodir,$R2Rstopdf_fn); # R-scape R2R PDF
  $R2Rstosvg = file($oodir,$R2Rstosvg_fn); # R-scape R2R SVG

#  open my $fh, ">", $out;
  my $rscape_options = "";
  if ($self->has_statistic){$rscape_options.=" --$self->statistic "}
  my $cmd = $rscape.$rscape_options.$self->ifile;
  my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmd, verbose => 0 );
  if( !$success ) {
    print STDERR "ERROR [$this_function] Call to $rscape unsuccessful\n";
    print STDERR "ERROR: $cmd\n";
    print STDERR "ERROR: this is what the command printed:\n";
    print join "", @$full_buf;
    croak $!;
  }
#  my $stdout_buffer = join "", @$stdout_buf;
#  my @rnaalifoldout = split /\n/, $stdout_buffer;
#  foreach my $line( @rnaalifoldout){
#    print $fh $line,"\n";
#  }
#  close($fh);

#  $self->_parse_rnaalifold($stdout_buffer);
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

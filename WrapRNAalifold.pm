# -*-CPerl-*-
# Last changed Time-stamp: <2017-01-22 19:54:57 mtw>

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
	       isa => 'MyDir3',
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

sub BUILD {
  my $self = shift;
  my $this_function = (caller(0))[3];
  confess "ERROR [$this_function] \$se+lf->alnfile not available"
    unless ($self->has_alnfile);
   $rnaalifold = can_run('RNAalfold') or
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
  my ($rnaalifold_outfilename,$rnaalifold_out);
  if ($self->has_alnfilebasename){$rnaalifold_outfilename = $self->alnfilebasename.".alifold.out"}
  elsif ($self->has_basename){$rnaalifold_outfilename = $self->bn.".alifold.out"}
  else{$rnaalifold_outfilename = "alifold.out"}
  $rnaalifold_out = file($oodir,$rnaalifold_outfilename);
  open my $fh, ">", $rnaalifold_out;
  my $rnaalifold_cmd = $rnaalifold." --aln --color -r --cfactor 0.6 --nfactor 0.5 -p ".$self->alnfile;
  my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $rnaalifold_cmd, verbose => 0 );
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
}

# parse RNAalifold output
sub _parse_rnaalifold {
  my ($self,$rnaz) = @_;
  my $this_function = (caller(0))[3];
  my @rnaz=split(/^/, $rnaz);
  my ($N,$identity,$columns,$decValue,$P,$z,$sci,$energy,$strand,
      $covariance,$combPerPair,$meanMFE,$consensusMFE,$consensusSeq,
      $consensusFold, $GCcontent, $ShannonEntropy);
  my @aln=();

  foreach my $i (0..$#rnaz){
    my $line=$rnaz[$i];
    $identity=$1 if ($line=~/Mean pairwise identity:\s*(-?\d+.\d+)/);
    $N=$1 if ($line=~/Sequences:\s*(\d+)/);
    if ($line=~/Reading direction:\s*(forward|reverse)/){
      $strand=($1 eq 'forward')?'+':'-';
    }
    $columns=$1 if ($line=~/Columns:\s*(\d+)/);
    $decValue=$1 if ($line=~/SVM decision value:\s*(-?\d+.\d+)/);
    $P=$1 if ($line=~/SVM RNA-class probability:\s*(-?\d+.\d+)/);
    $z=$1 if ($line=~/Mean z-score:\s*(-?\d+.\d+)/);
    $sci=$1 if ($line=~/Structure conservation index:\s*(-?\d+.\d+)/);
    $energy=$1 if ($line=~/Energy contribution:\s*(-?\d+.\d+)/);
    $covariance=$1 if ($line=~/Covariance contribution:\s*(-?\d+.\d+)/);
    $combPerPair=$1 if ($line=~/Combinations\/Pair:\s*(-?\d+.\d+)/);
    $consensusMFE=$1 if ($line=~/Consensus MFE:\s*(-?\d+.\d+)/);
    $meanMFE=$1 if ($line=~/Mean single sequence MFE:\s*(-?\d+.\d+)/);
    $GCcontent=$1 if ($line=~/G\+C content:\s(\d+.\d+)/);
    $ShannonEntropy=$1 if ($line=~/Shannon entropy:\s*(\d+.\d+)/);

    if ($line=~/^>/){
      chomp($rnaz[$i+1]);
      chomp($rnaz[$i+2]);
      if ($line=~/^>consensus/){
	$consensusSeq=$rnaz[$i+1];
	$consensusFold=substr($rnaz[$i+2],0,length($rnaz[$i+1]));
	last;
      } else {
	if ($line=~/>(.*?) (\d+) (\d+) (\+|\-) (\d+)/){
	  push @aln, {name=>$1,
		      start=>$2,
		      end=>$2+$3,
		      strand=>$4,
		      fullLength=>$5,
		      seq=>$rnaz[$i+1],
		      fold=>substr($rnaz[$i+2],0,length($rnaz[$i+1]))};
	  $i+=2;
	} elsif ($line=~/^(.*)\/(\d+)-(\d+)$/){
	  push @aln, {name=>$1,
		      start=>$2,
		      end=>$3,
		      strand=>$strand,
		      fullLength=>'',
		      seq=>$rnaz[$i+1],
		      fold=>substr($rnaz[$i+2],0,length($rnaz[$i+1]))};
	  $i+=2;
	}
      }
    }
  }

  $self->P($P);
  $self->z($z);
  $self->sci($sci);
}

1;

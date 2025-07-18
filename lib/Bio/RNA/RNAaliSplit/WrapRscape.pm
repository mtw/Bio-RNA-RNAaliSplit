# -*-CPerl-*-
# Last changed Time-stamp: <2025-07-15 16:17:33 mtw>
# place of birth: somewhere over Newfoundland

# Bio::RNA::RNAaliSplit::WrapRscape.pm: A versatile object-oriented
# wrapper for R-scape
#
# R-scape output files are created in a temporary directory below the
# object's output directory and then moved to their final destination.
# This avoids collisions when multiple instances run in the same
# directory.
#
# Requires R-scape v2 or above or above available to the Perl
# interpreter.

package Bio::RNA::RNAaliSplit::WrapRscape;

use version; our $VERSION = qv('0.12');
use Carp;
use Data::Dumper;
use Moose;
use Path::Class;
use IPC::Cmd qw(can_run run);
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use File::pushd;
use Cwd;

my ($rscape,$oodir);
my $exe = "R-scape";

has 'basename' => (
		   is => 'rw',
		   isa => 'Str',
		   predicate => 'has_basename',
		   documentation => q(Set this to override output basename),
		  );

has 'statistic' => (
		    is => 'rw',
		    isa => 'Str',
		    predicate => 'has_statistic',
		    default => 'GTp',
		    documentation => q(Covariation statistic),
		   );
has 'cseq' => (
		    is => 'rw',
		    isa => 'Int',
		    predicate => 'has_cseq',
		    documentation => q(Number of sequences counted from stk file),
		    init_arg => undef,
	      );

has 'nseq' => (
		    is => 'rw',
		    isa => 'Int',
		    predicate => 'has_nseq',
		    documentation => q(Number of sequences returned by R-scape),
		    init_arg => undef,
	      );

has 'alen' => (
		    is => 'rw',
		    isa => 'Int',
		    predicate => 'has_alen',
		    documentation => q(Alignment length),
		    init_arg => undef,
	      );

has 'nbpairs' => (
		    is => 'rw',
		    isa => 'Int',
		    predicate => 'has_nbpairs',
		    documentation => q(Number of base pairs),
		    init_arg => undef,
		 );

has 'evalue' => (
		    is => 'rw',
		    isa => 'Num',
		    predicate => 'has_evalue',
		    documentation => q(E-value threshold),
		    init_arg => undef,
		);

has 'FP' => (
		    is => 'rw',
		    isa => 'Int',
		    predicate => 'has_FP',
		    documentation => q(Covarying non base pairs (FP)),
		    init_arg => undef,
	    );

has 'TP' => (
	     is => 'rw',
	     isa => 'Int',
	     predicate => 'has_TP',
	     documentation => q(Covarying base pairs (TP)),
	     init_arg => undef,
	    );

has 'T' => (
	     is => 'rw',
	     isa => 'Int',
	     predicate => 'has_T',
	     documentation => q(Number of base pairs (T)),
	     init_arg => undef,
	   );

has 'F' => (
	    is => 'rw',
	    isa => 'Int',
	    predicate => 'has_F',
	    documentation => q(Total number of covarying base pairs (F)),
	    init_arg => undef,
	   );

has 'Sen' => (
	      is => 'rw',
	      isa => 'Num',
	      predicate => 'has_Sen',
	      documentation => q(Sensitivity (TP/T)),
	      init_arg => undef,
	     );

has 'PPV' => (
	      is => 'rw',
	      isa => 'Num',
	      predicate => 'has_PPV',
	      documentation => q(Positive predictive value (TP/F)),
	      init_arg => undef,
	     );

has 'Fmeasure' => (
		   is => 'rw',
		   isa => 'Num',
		   predicate => 'has_Fmeasure',
		   documentation => q(F-measure (2*Sen*PPV(Sen+PPV))),
		   init_arg => undef,
	     );

has 'nofigures' => (
		    is => 'rw',
		    isa => 'Int',
		    predicate => 'has_nofigures',
		    documentation => q(Turn off all image procudtion by R-scape),
		   );

has 'naive' => (
		    is => 'rw',
		    isa => 'Int',
		    predicate => 'has_naive',
		    documentation => q(Turn on R-scape naive mode),
		   );

has 'sigBP' => ( # significantly covarying base pairs
		is => 'rw',
		isa => 'ArrayRef',
		default => sub { [] },
		traits => ['Array'],
		predicate => 'has_data',
		handles => {
			    all    => 'elements',
			    count  => 'count',
			    add    => 'push',
			    pop    => 'pop',
			   },
	       );

has 'status' => (
		 is => 'rw',
		 isa => 'Int',
		 predicate => 'has_status',
		 documentation => q(R-scape program status), # 0:OK ; >0:NOTOK
		);


with 'FileDirUtil';

sub BUILD {
  my $self = shift;
  my $this_function = (caller(0))[3];
  confess "ERROR [$this_function] \$self->ifile not available"
    unless ($self->has_ifile);

  # locate R-scape executable
  $rscape = can_run($exe) or
    croak "ERROR [$this_function] $exe not found";

  # set default output directory if none provided
  unless($self->has_odir){
    my $odir_name = "as";
    $self->odir( [$self->ifile->dir,$odir_name] );
  }

  # enforce naive mode for specific statistics
  if ($self->has_statistic
      && $self->statistic =~ /^(?:RAF|RAFS|RAFp|RAFa|RAFSp|RAFSa)$/) {
    $self->naive(1);
  }

  # prepare output subdirectory
  $oodir = $self->odir->subdir("rscape");
  my @created = make_path($oodir, {error => \my $err});
  confess "ERROR [$this_function] could not create output directory $$oodir"
    if (@$err);

  $self->set_ifilebn;
  $self->_count_seq();
  if ($self->cseq > 1){ $self->run_rscape() }
}

sub _count_seq {
  my $self = shift;
  my $this_function = (caller(0))[3];
  my $count = 0;
  open my $stk, "<", $self->ifile
    or croak "ERROR [$this_function] Cannot open Stockholm file $self->ifile for reading";
  while (<$stk>){
    next if (/^#/);
    next if (/^\/\//);
    $count++;
  }
  close ($stk);
  $self->cseq($count);
}

sub run_rscape {
  my $self = shift;
  my $this_function = (caller(0))[3];
  my ($out_fn,$sout_fn,$out,$sout,$sum);
  my ($rscape_out_bn,$rscape_out,$rscape_sout,$rscape_sum);
  my ($tmpdir, $dir, $cwd);
  my $tag = "";
  my $tmpprefix = '';
  my $ifile_abs_path = $self->ifile->absolute->stringify;
  if ($self->has_statistic){$tag = ".".$self->statistic};

  # work in a temporary directory inside odir to avoid file name clashes
  if ($self->has_basename){
    $tmpprefix = $self->bn."_";
  }
  elsif ($self->has_ifilebn){
    $tmpprefix = $self->ifilebn."_";
  }
  else{
    $tmpprefix = "XXXX_";
  }
  $tmpprefix .= "XXXXX";
  $tmpdir = tempdir( $tmpprefix, DIR => $oodir );
  my @created = make_path($tmpdir, {error => \my $err});
  confess "ERROR [$this_function] could not create output directory $oodir"
    if (@$err);

  $dir = pushd($tmpdir);
  $cwd = getcwd();

  if ($self->has_basename){
    $out_fn  = $self->basename.$tag."."."rscape.out";
    $sout_fn = $self->basename.$tag."."."rscape.sorted.out";
  }
  elsif ($self->has_ifilebn){
    $out_fn  = $self->ifilebn.$tag."."."rscape.out";
    $sout_fn = $self->ifilebn.$tag."."."rscape.sorted.out";
  }
  else{
    $out_fn  = $tag."rscape.out";
    $sout_fn = $tag."rscape.sorted.out";
  }
  $out  = file($out_fn);  # R-scape stdout
  $sout = file($sout_fn); # R-scape sorted stdout

  $rscape_out_bn = "rscape.out";
  $rscape_out  = file($rscape_out_bn); # produced inside tmpdir
  $rscape_sout = file($rscape_out_bn.".sorted");

  my $rscape_options = " --outname $rscape_out_bn --rna --outdir . ";
  if ($self->has_nofigures && $self->nofigures == 1){$rscape_options.=" --nofigures "};
  if ($self->has_naive && $self->naive == 1){$rscape_options.=" --naive "};
  if ($self->has_statistic){$rscape_options.=" --".$self->statistic." "  }
  my $cmd = $rscape.$rscape_options.$ifile_abs_path;

  my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmd, verbose => 0 );
  if( !$success ) {
    print STDERR "ERROR [$this_function] Call to $rscape unsuccessful\n";
    print STDERR "ERROR: $cmd\n";
    print STDERR "ERROR: this is what the command printed:\n";
    print join "", @$full_buf;
    croak $!;
  }

  open my $fh, ">", $out
    or croak "ERROR: cannot open '$out' for writing: $!\n";
  my $stdout_buffer = join "", @$stdout_buf;
  my @out_buffer = split /\n/, $stdout_buffer;
  foreach my $line( @out_buffer){
    print $fh $line,"\n";
  }
  close($fh);

  $self->_parse_rscape($stdout_buffer);

  rename $rscape_out, $out;
  rename $rscape_sout, $sout;
}

# parse R-scape output
# R-scape output is expected to look like this:
#
#  MSA myaln nseq 9 (9) alen 45 (45) avgid 79.51 (79.51) nbpairs 11 (11)
#  contacts  11 (11 bpairs 11 wc bpairs)
#  maxD      8.00
#  mind      1
#  Method Target_E-val [cov_min,conv_max] [FP | TP True Found | Sen PPV F]
#  RAFS    0.05           [-1.00,0.83]    [0 | 3 11 3 | 27.27 100.00 42.86]
#        left_pos       right_pos        score   E-value
# ------------------------------------------------------------
#*	        20	        36	0.83	0.000328218
#*	        21	        35	0.83	0.000328218
#*	        22	        34	0.33	0.0497359

sub _parse_rscape {
  my ($self,$out) = @_;
  my $this_function = (caller(0))[3];
  my @buffer=split(/^/, $out);
  my $parse1 = 0;
  my $parse2 = 0;
  my $nosbp = 0;

  foreach my $i (0..$#buffer){
    chomp($buffer[$i]);

    if ($buffer[$i] =~ m/^#\s+MSA\s+([a-zA-Z0-9_.|]+)\s+nseq\s+(\d+)\s+\(\d+\)\s+alen\s+(\d+)\s+\(\d+\)\s+avgid\s+(\d+\.\d+)\s+\(\d+\.\d+\)\s+nbpairs\s+(\d+)\s+\(\d+\)/g){
      $self->nseq($2);
      $self->alen($3);
      $self->nbpairs($5);
      $parse1 = 1;
      next;
    }
    if ($buffer[$i] =~ m/^#\s+([a-zA-Z0-9]+)\s+(\d+\.\d+)\s+\[(\-?\d+\.\d+),(-?\d+\.\d+)\]\s+\[(\d+)\s+\|\s+(\d+)\s+(\d+)\s+(\d+)\s+\|\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\]/g){
      $self->evalue($2);
      $self->FP($5);
      $self->TP($6);
      $self->T($7);
      $self->F($8);
      $self->Sen($9);
      $self->PPV($10);
      $self->Fmeasure($11);
      $parse2 = 1;
      next;
    }
    if ($buffer[$i] =~ m/^\*\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)$/){
      my %bp = (i => $1,
		j => $2,
		score => $3,
		evalue => $4);
      push @{$self->sigBP}, \%bp;
    }
    if ($buffer[$i] =~ m/^no significant pairs$/){
      $nosbp=1;
    }
  }

  #carp "INFO: [$this_function] parse1:".eval($parse1);
  #carp "INFO: [$this_function] parse2:".eval($parse2);

  if ($nosbp == 1){
    $self->status(1); # no significant basepairs
    $self->_check_attributes();
    return;
  }
  if ($parse1 == 1 && $parse2 == 1){
    $self->_check_attributes();
    $self->status(0); # all OK
    return;
  }
  elsif ($parse1 == 1 && $parse2 == 0){
    $self->status(2); # covariation scores are almost constant, no further analysis
    return;
  }
  else{
    croak "ERROR: [$this_function] ambiguous status when parsing rscape output file. parse1:".eval($parse1)." parse2:".eval($parse2)." This shouldn't happen ...";
    $self->status(3);
    return;
  }
}

sub _check_attributes {
  my ($self,$out) = @_;
  my $this_function = (caller(0))[3];

  $self->has_cseq     ? 1 : croak "ERROR [$this_function] \$self->cseq not set ". $self->ifile;
  $self->has_nseq     ? 1 : croak "ERROR [$this_function] \$self->nseq not set ". $self->ifile;
  $self->has_alen     ? 1 : croak "ERROR [$this_function] \$self->alen not set". $self->ifile;
  $self->has_nbpairs  ? 1 : croak "ERROR [$this_function] \$self->nbpairs not set". $self->ifile;
  $self->has_evalue   ? 1 : croak "ERROR [$this_function] \$self->evalue not set". $self->ifile;
  $self->has_FP       ? 1 : croak "ERROR [$this_function] \$self->FP not set". $self->ifile;
  $self->has_TP       ? 1 : croak "ERROR [$this_function] \$self->TP not set". $self->ifile;
  $self->has_T        ? 1 : croak "ERROR [$this_function] \$self->T not set". $self->ifile;
  $self->has_F        ? 1 : croak "ERROR [$this_function] \$self->F not set". $self->ifile;
  $self->has_Sen      ? 1 : croak "ERROR [$this_function] \$self->Sen not set". $self->ifile;
  $self->has_PPV      ? 1 : croak "ERROR [$this_function] \$self->PPV not set". $self->ifile;
  $self->has_Fmeasure ? 1 : croak "ERROR [$this_function] \$self->Fmeasure not set". $self->ifile;
}

1;

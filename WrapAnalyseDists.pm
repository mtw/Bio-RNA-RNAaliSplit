# -*-CPerl-*-
# Last changed Time-stamp: <2017-01-24 19:36:17 mtw>

# WrapAnalyseDists.pm: Wrapper for computing split decomposition
#
# requires AnalyseDists executable from the ViennaRNA poackage
# available to the Perl interpreter

package WrapAnalyseDists;

use version; our $VERSION = qv('0.01');
use Carp;
use Data::Dumper;
use Moose;
use Moose::Util::TypeConstraints;
use Path::Class::File;
use Path::Class::Dir;
use Path::Class;
use IPC::Cmd qw(can_run run);
use Array::Set qw(set_diff);

my ($analysedists,$oodir);

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

has 'splits' => (
		 is => 'rw',
		 isa => 'ArrayRef',
		 default => sub { [] },
		 predicate => 'has_splits',
		 traits => ['Array'],
		 handles => {
			     allsplits => 'elements',
			     count     => 'count',
			     add       => 'push',
			     pop       => 'pop',
			    },
		);

has 'nr_splits' => (
		    is => 'rw',
		    isa => 'Num',
		    predicate => 'has_nr_splits',
		   );

has 'dim' => (
	      is => 'rw',
	      isa => 'Num',
	      predicate => 'has_dim',
	     );

with 'FileDir';

sub BUILD {
  my $self = shift;
  my $this_function = (caller(0))[3];
  confess "ERROR [$this_function] \$self->ifile not available"
    unless ($self->has_ifile);
   $analysedists = can_run('AnalyseDists') or
    croak "ERROR [$this_function] AnalyseDists not found";
  unless($self->has_odir){
    unless($self->has_odirname){
      self->odirname("as");
    }
    $self->odir( [$self->ifile->dir,$self->odirname] );
    mkdir($self->odir);
  }
  $oodir = $self->odir->subdir("analysedists");
  mkdir($oodir);
  $self->dim( $self->_get_dim() );

  # do computation
  $self->_NeighborJoining();
  $self->SplitDecomposition();
  $self->nr_splits($self->count);
}

sub _NeighborJoining {
  # TODO  warn if negative branch lengths occur
  my $self = shift;
  my $this_function = (caller(0))[3];
  my ($nj_outfilename,$nj_treefilename,$nj_out,$nj_tree);

  if ($self->has_basename){
    $nj_outfilename = $self->basename.".nj.out";
    $nj_treefilename = $self->basename.".nj.ps";
  }
  else{
    $nj_outfilename = "nj.out";
    $nj_treefilename = "nj.ps";
  }
  $nj_out = file($oodir,$nj_outfilename);
  $nj_tree = file($oodir,$nj_treefilename);
  open my $fh, ">", $nj_out;

  my $ad_cmd = $analysedists." -Xn < ".$self->ifile;
  my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $ad_cmd, verbose => 0 );
  if( !$success ) {
    print STDERR "ERROR [$this_function] Call to $analysedists unsuccessful\n";
    print STDERR "ERROR: this is what the command printed:\n";
    print join "", @$full_buf;
    croak $!;
  }
  my $stdout_buffer = join "",@$stdout_buf;
  my @out = split /\n/, $stdout_buffer;
  foreach my $line( @out){print $fh $line,"\n"}
  close($fh);
  rename "nj.ps", $nj_tree;
  $self->_parse_nj($stdout_buffer);
}

# parse the output of AnalyseDists -Xn
# populate array of hashes, each holding two sets of nodes corresponding to splits
sub _parse_nj {
  my ($self,$nj) = @_;
  my $this_function = (caller(0))[3];
  my %data = ();
  my $num;
  my $count = 1;
  my @lines =  split /\n/,$nj;
  foreach my $line (@lines){
    my @s1 = ();
    my @set1 = ();
    my @set2 = ();
    next if ($line =~ m/^>\s+\D/);
    if ($line =~ m/^>\s+(\d+)/){$num = $1;next}
    last if ($count++ >= $num);
    my @all = (1..$num);
    print " #### $line\n";
    croak "ERROR [$this_function] Cannot parse neighbor joining graph line\n$line\n"
      unless ($line =~ m/^\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)/g);
    my $i = $1;
    my $j = $2;
    print "$i -> $j\n";
    push @{$data{$i}}, $j;
    if (exists $data{$j}){push @{$data{$i}}, @{$data{$j}} };
    print Dumper(\%data);
    # populate set1
    push @s1, $i;
    push @s1, @{$data{$i}};
    @set1 =  sort {$a <=> $b} @s1;
    my @diff =  set_diff(\@all, \@set1);
    @set2 =  sort {$a <=> $b} @{$diff[0]};
    print Dumper(\@set1);
    print Dumper(\@set2);
    print "+++++++++++++++++++++++++++++++++++\n";
  }
}

sub SplitDecomposition {
  my $self = shift;
  my $this_function = (caller(0))[3];
  my ($sd_outfilename,$sd_out);
  if ($self->has_basename){$sd_outfilename = $self->basename.".sd.out"}
  else{$sd_outfilename = "sd.out"}
  $sd_out = file($oodir,$sd_outfilename);
  open my $fh, ">", $sd_out;

  my $sd_cmd = $analysedists." -Xs < ".$self->ifile;
  my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $sd_cmd, verbose => 0 );
  if( !$success ) {
    print STDERR "ERROR [$this_function] Call to $analysedists unsuccessful\n";
    print STDERR "ERROR: this is what the command printed:\n";
    print join "", @$full_buf;
    croak $!;
  }
  my @sdout = split /\n/, $$stdout_buf[0];
  foreach my $line( @sdout){
    print $fh $line,"\n";
  }
  close($fh);
  $self->_parse_sd($$stdout_buf[0]); # parse split graph data
}

# parse the output of AnalyseDists -Xs
# populate array of hashes, each holding two sets of nodes corresponding to splits
sub _parse_sd {
  my ($self,$sd) = @_;
  my $this_function = (caller(0))[3];
  my $num;
  my @lines =  split /\n/,$sd;
  foreach my $line (@lines){
    next if ($line =~ m/^>\s+\D/);
    if ($line =~ m/^>\s+(\d+)/){$num = $1;next}
    last if ($line =~ m/^\s+\d+\.\d+\s+\:\s+\{\s+\[Split prime fraction\]\s+\}/g );
 #   print "$line\n";
    croak "ERROR [$this_function] Cannot parse split graph line\n$line\n"
      unless ($line =~ m/^\s+\d+\s+\d+\.\d+\s+:\s+\{\s+([\d+\s+]+)\|/g);
    my @foo = split /\s+/, $1; # set 1
    my @moo = (1 .. $self->dim);
    my @bar = (); # set 2
    foreach my $i (@moo){
      push (@bar, $i) unless ( grep {$i == $_}@foo );
    }
    $self->add( {S1=>\@foo,S2=>\@bar } );
  }
  croak "ERROR [$this_function] expected $num splits, but parsed ".$self->count
    unless ($self->count == $num);
}

sub _get_dim {
  my $self = shift;
  my $this_function = (caller(0))[3];
  my $dim = -1 ;
  open my $fh, "<", $self->ifile or die $!;
  while(<$fh>){
    if (m/^>\s+X\s+(\d+)/){$dim = $1;last;}
  }
  croak "ERROR [$this_function] could not parse dimension from input matrix"
    if ($dim == -1);
  close($fh);
  return $dim;
}

1;


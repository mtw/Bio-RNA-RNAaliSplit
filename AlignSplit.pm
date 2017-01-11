# -*-CPerl-*-
# Last changed Time-stamp: <2017-01-11 19:58:32 mtw>
package AlignSplit;

use version; our $VERSION = qv('0.01');
use Carp;
use Data::Dumper;
use Moose;
use Moose::Util::TypeConstraints;
use Path::Class::File;
use Path::Class::Dir;
use Path::Class;
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

subtype 'MyDir' => as class_type('Path::Class::Dir');

coerce 'MyDir'
  => from 'ArrayRef'
  => via { Path::Class::Dir->new( @{ $_ } ) };

has 'infile' => (
	       is => 'ro',
	       isa => 'MyFile',
	       predicate => 'has_infile',
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

has 'next_aln' => (
		   is => 'rw',
		   isa => 'Bio::SimpleAlign',
		   predicate => 'has_next_aln',
		   init_arg => undef, # make this unsettable via constructor
		  );

has 'odirname' => (
		   is => 'ro',
		   default => 'as',
		   predicate => 'has_odirname',
		  );

has 'odir' => (
	       is => 'rw',
	       isa => 'MyDir',
	       predicate => 'has_odir',
	       coerce => 1,
	      );

has 'dump' => (
	       is => 'rw',
	       isa => 'Num',
	       predicate => 'has_dump_flag',
	       );

has 'sci' => (
	      is => 'rw',
	      isa => 'Num',
	      predicate => 'has_sci',
	     );

has 'consensus_struc' => (
			  is => 'rw',
			  isa => 'Str',
			  predicate => 'has_consensus_struc',
			 );

has 'consensus_mfe' => (
			is => 'rw',
			isa => 'Num',
			predicate => 'has_consensus_mfe',
		       );
has 'consensus_energy' => (
			   is => 'rw',
			   isa => 'Num',
			   predicate => 'has_consensus_energy',
			  );

has 'consensus_covar_terms' => (
				is => 'rw',
				isa => 'Num',
				predicate => 'has_consensus_covar_terms',
			       );

sub BUILD {
    my $self = shift;
    my $this_function = (caller(0))[3];
    my $odir;
    confess "ERROR [$this_function] \$self->infile not available"
	unless ($self->has_infile);
    $self->alignment({
		      -file => $self->infile,
		      -format => $self->format});
    $self->next_aln($self->alignment->next_aln);
    $self->odir( [$self->infile->dir,$self->odirname ] );
    mkdir($self->odir);
    if ($self->has_dump_flag){ # dump input to odir/input
      my $iodir = $self->odir->subdir('input');
      mkdir($iodir);
      _dump_infile_as_aln($self,$iodir);
   #   dump_infile_as_seq($iodir);
      print "HAS DUMP FLAG\n";
    }
    $self->_alifold();
  }

sub dump_subalignment {
  # TODO receive array indicating which sequences to dump into subalignment
  my ($self,$alipathsegment,$alinr) = @_;
  my $this_function = (caller(0))[3];

  # create subalignment output path and file
  my $ids = join "_", @$alinr;
  unless (defined($alipathsegment)){$alipathsegment = "tmp"}
  my $oodir = $self->odir->subdir($alipathsegment);
  mkdir($oodir);
  my $oalifile = file($oodir,$ids.".aln");
  $oalifile->touch;

  my $outali = Bio::AlignIO->new(-file   => ">$oalifile",
#				 -fh => \*STDOUT, 
				 -format => "ClustalW",
				 -flush  => 0 ); # go as fast as we can!
  #  print Dumper($outali);

  my $aln2 = $self->next_aln->select_noncont(@$alinr);
#  print Dumper($aln2);
  $outali->write_aln( $aln2 );
  return ( $oalifile );
}


sub _alifold {
  my $self = shift;
  my $this_function = (caller(0))[3];
  my $sci = 0.0;
  my $f;
  my $alifold = can_run('RNAalifold') or
    croak " ERROR [$this_function] RNAalifold not found";
 # carp "[$this_function]: $alifold";

  $f = $self->infile->stringify;
  my $cmd = "$alifold --sci $f";

  my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmd, verbose => 0 );
  if( !$success ) {
    print STDERR "ERROR [$this_function] Call to $alifold unsuccessful\n";
    print STDERR "ERROR: this is what the command printed:\n";
    print join "", @$full_buf;
    croak $!;
  }
  my @aliout = split /\n/, $$stdout_buf[0];
  croak "ERROR [$this_function] cannot parse RNAalifold output"
    unless ($aliout[1] =~ m/([\(\)\.]+)\s+\((-?\d+\.\d+)\s+=\s+(-?\d+\.\d+)\s+\+\s+(-?\d+\.\d+)\)\s+\[sci\s+=\s+(\d+\.\d+)\]/);
  $self->consensus_struc($1);
  $self->consensus_mfe($2);
  $self->consensus_energy($3);
  $self->consensus_covar_terms($4);
  $self->sci($5);
}


sub _dump_infile_as_aln {
  my ($self,$dir) = @_;
  my $this_function = (caller(0))[3];
  print ref($self);
  #  print Dumper($self->next_aln);
  #foreach my $seq ($self->next_aln->each_seq) {
  #  print ">> $seq->seq\n";
  #}
}

no Moose;

1;

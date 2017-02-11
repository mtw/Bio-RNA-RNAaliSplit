# -*-CPerl-*-
# Last changed Time-stamp: <2017-02-12 00:27:52 mtw>

# FileDir.pm: Moose Role for basic file IO

package FileDir;

use version; our $VERSION = qv('0.01');
use Moose::Util::TypeConstraints;
use Moose::Role;
use Path::Class::File;
use Path::Class::Dir;

subtype 'MyFile' => as class_type('Path::Class::File');

coerce 'MyFile'
  => from 'Str'
  => via { Path::Class::File->new($_) };

subtype 'MyDir' => as class_type('Path::Class::Dir');

coerce 'MyDir'
  => from 'ArrayRef'
  => via { Path::Class::Dir->new( @{ $_ } ) };

has 'ifile' => (
		is => 'ro',
		isa => 'MyFile',
		predicate => 'has_ifile',
		coerce => 1,
	      );

has 'odir' => (
	       is => 'rw',
	       isa => 'MyDir',
	       predicate => 'has_odir',
	       coerce => 1,
	      );

has 'odirn' => ( # custom output dir name
		is => 'rw',
		isa => 'Path::Class::Dir',
		predicate => 'has_odirn',
	       );


no Moose;

1;

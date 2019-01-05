# -*-CPerl-*-
# Last changed Time-stamp: <2019-01-05 21:28:20 mtw>

package Bio::RNA::RNAaliSplit::Roles;

use version; our $VERSION = qv('0.08');
use Moose::Util::TypeConstraints;
use Moose::Role;
use Path::Class::Dir;
use namespace::autoclean;

has 'dirnam' => ( # custom output dir name
		  is => 'rw',
		  isa => 'Path::Class::Dir',
		  predicate => 'has_dirnam',
		 );

1;

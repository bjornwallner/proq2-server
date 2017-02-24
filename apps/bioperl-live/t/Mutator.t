# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
my $error;
use vars qw($NUMTESTS);
BEGIN { 
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    $error=0;
    use Test;
    $NUMTESTS=12;
    plan tests => $NUMTESTS;
    eval { require 'IO/String.pm' };
    if( $@ ) {
	print STDERR "IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping tests.\n";
	for( 1..$NUMTESTS ) {
	    skip(1,"IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping tests");
	}
	$error = 1; 
    }

}

if( $error ==  1 ) {
    exit(0);
}

require Bio::LiveSeq::Mutator;
require Bio::LiveSeq::IO::BioPerl;
require Bio::LiveSeq::Gene;
require Bio::Root::IO;


$a = Bio::LiveSeq::Mutator->new();
ok $a;

ok $a->numbering, 'coding';
ok $a->numbering('coding 1');
ok $a->numbering, 'coding 1';

require Bio::LiveSeq::Mutation;
my $mt = new Bio::LiveSeq::Mutation;
ok $mt->seq('g');
$mt->pos(100);
ok ($a->add_Mutation($mt));
my @each = $a->each_Mutation;
ok( (scalar @each), 1 );
my $mt_b = pop @each;
ok($mt_b->seq, 'g');
my $filename=Bio::Root::IO->catfile("t","data","ar.embl");
my $loader=Bio::LiveSeq::IO::BioPerl->load('-file' => "$filename");
my $gene_name='AR'; # was G6PD

my $gene=$loader->gene2liveseq('-gene_name' => $gene_name, 
			       '-getswissprotinfo' => 0);
ok($gene);
ok $a->gene($gene);

my $results = $a->change_gene();

ok($results);
use Bio::Variation::IO;
if ($results) {    
    my $out = Bio::Variation::IO->new( '-format' => 'flat');
    ok($out->write($results));    
}

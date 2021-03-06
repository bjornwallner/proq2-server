# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id: RootI.t,v 1.3 2001/01/25 22:13:40 jason Exp $

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;    
    plan tests => 7;
}

use Bio::Root::RootI;

my $obj = new Bio::Root::RootI();

ok defined($obj) && $obj->isa('Bio::Root::RootI');

eval { $obj->throw('Testing throw') };
ok $@ =~ /Testing throw/;# 'throw failed';

# doesn't work in perl 5.00405
#my $val;
#eval {
#    my ($tfh,$tfile) = $obj->tempfile();
#    local * STDERR = $tfh;
#    $obj->warn('Testing warn');
#    close $tfh;    
#    open(IN, $tfile) or die("cannot open $tfile");    
#    $val = join("", <IN>) ;
#    close IN;
#    unlink $tfile;
#};
#ok $val =~ /Testing warn/;
#'verbose(0) warn did not work properly' . $val;

$obj->verbose(-1);
eval { $obj->throw('Testing throw') };
ok $@=~ /Testing throw/;# 'verbose(-1) throw did not work properly' . $@;

eval { $obj->warn('Testing warn') };
ok !$@;

$obj->verbose(1);
eval { $obj->throw('Testing throw') };
ok $@ =~ /Testing throw/;# 'verbose(1) throw did not work properly' . $@;

# doesn't work in perl 5.00405
#undef $val;
#eval {
#    my ($tfh,$tfile) = $obj->tempfile();
#    local * STDERR = $tfh;
#    $obj->warn('Testing warn');
#    close $tfh;
#    open(IN, $tfile) or die("cannot open $tfile");    
#    $val = join("", <IN>);
#    close IN;
#    unlink $tfile;
#};
#ok $val =~ /Testing warn/;# 'verbose(1) warn did not work properly' . $val;

my @stack = $obj->stack_trace();
ok scalar @stack, 2;

my $verbobj = new Bio::Root::RootI(-verbose=>1,-strict=>1);
ok $verbobj->verbose(), 1;

1;



#! /usr/bin/perl

use utf8;
use strict;
use warnings;
use line_svg;

my $sub_cmd = shift;
my %cmds = map{$_,1} qw(line_svg);

unless (defined $sub_cmd) { die help_text(); };
unless (exists $cmds{$sub_cmd}) {
	warn ' Please give valid sub command ! ', "\n";
	die help_text();
}

SWITCH:{
	$sub_cmd eq 'line_svg' && do { 
		my $self = line_svg->new();
		$self->process(\@ARGV);
		last SWITCH;
		};
	$sub_cmd eq 'help'   && do { die help_text(); last SWITCH; };
}

sub help_text {
	return <<HELP
	Usage:   graph_main.pl <command> [options]
		line_svg    the graph like is lines
		xxxx        thexxxx
		xxxx        thexxx
HELP
}

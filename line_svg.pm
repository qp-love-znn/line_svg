#!/usr/bin/perl -w

package line_svg;
require Exporter;
use strict;
use warnings;
use SVG;
use Getopt::Long;

our @ISA = qw(Exporter);
our @EXPORT = qw(circ_gene);
our @version = 1.0;

sub new{
	my $class = shift;
	my $self = {};
	bless $self;
#	$self->{_CIRC_POSIT} = undef;
#	$self->{_GNENE_POSIT} = undef;
	$self->{_CONF} = undef;
	$self->{_TYPE} = undef;
	$self->{_OUTPUT} = undef;
	return $self;
}

sub process {
	my $self = shift;
	my ( $help, $options );
	unless( @ARGV ) { die $self->help_text(); };
	$options = GetOptions (
		'conf=s' => \$self->{_CONF},
		'type=s' => \$self->{_TYPE},
		'output=s' => \$self->{_OUTPUT},
		'help' => \$help,
		);
	if ( $help ) { print STDERR help_text(); exit 0; };
	# Check on all the input data
	print STDERR "input file not found or is empty: $self->{_CONF}\n" unless( -s $self->{_CONF} );
	print STDERR "Output directory not found: $self->{_OUTPUT}\n" unless($self->{_OUTPUT});
	print STDERR "please give the graph type  $self->{_TYPE}\n" unless($self->{_TYPE});
	if($self->{_TYPE} eq "circ_gene"){
		$self->circ_gene($self->{_CONF},$self->{_OUTPUT});
	}
}

##get the graph for GHR160145-sup8##
sub circ_gene{
	my $self = shift;
	my ($conf,$out) = @_;
	open IN,$conf or die $!;
	print "read the circ_gene.conf\n";
	my %conf_hash;
	while(<IN>){
		chomp;
		my @t = split /\s+=\s+/,$_;
		$conf_hash{$t[0]} = $t[1];
	}
	close IN;
	open CIR,$conf_hash{"cir_posit"} or die $!;
	my (%cir_s,%cir_e,%genes_s,%genes_e);
	my $cir_n = 0;#the number of circ;
	print "read the circ_posit\n";
	while(<CIR>){
		chomp;
		my @t = split /\t/,$_;
		$cir_n += 1;
		$cir_s{$t[0]} = $t[1];
		$cir_e{$t[0]} = $t[2];
	}
	close CIR;
	my @col = ("#191970","#43CD80","#FFC125","#CD7054","#66CD00","#FFA500","#9370DB","#9400D3","#FF4500","#FF6347","#F4A460");
	my @g = ("misc_RNA","NS1/1C","NS2/1B","NP","P","M","SH","G","F0","22K/M2","L");
	my %gcol;
	for my $i (0..$#g){
		$gcol{$g[$i]} = $col[$i];
	}
	open IN1,$conf_hash{"gene_posit"} or die $!;
	my $n = 0;
	my $genes_n = 0;#the number of gene;
	while(<IN1>){
		chomp;
		my @t = split /\t/,$_;
		$genes_n += 1;
		$genes_s{"$t[1]\t$t[2]"} = $t[0];
	}
	close IN1;
	my $svg = SVG->new(width => $conf_hash{"width"},height => $conf_hash{"height"});
	my $forward_y = 30;#上顶端留10
	my $forward_x = 100;#最左端留10
	my $inter_y = 0;#每条环状为一条线，间隔为1
	for my $i (0..19){
		my $x = $i*50;
		$svg->line(x1=>$x+100,y1=>22,,x2=>100+$x,y2=>27,stroke=>"black",strokewidth=>1);
	}
	$svg->line(x1=>1100,y1=>22,,x2=>1100,y2=>27,stroke=>"black",strokewidth=>1);
	$svg->line(x1=>100,y1=>27,x2=>100,y2=>27,x2=>100+1000,stroke=>"black",strokewidth=>1);
	for my $c (keys %cir_s){
		my $x1 = (1000*$cir_s{$c})/$conf_hash{"chr_length"};#15222为chr的长度，1000为设定的长度
		my $x2 = (1000*$cir_e{$c})/$conf_hash{"chr_length"};
		$svg->line(x1=>$x1+$forward_x,y1=>$forward_y+$inter_y,x2=>$x2+$forward_x,y2=>$forward_y+$inter_y,stroke=>"green",strokewidth=>0.5);
		$inter_y += 0.8;#每条环状为一条线，间隔为1
	}
	my $gene_y = $forward_y + $inter_y + 5;
	my $inter_gene_y = 0;
	for my $po (keys %genes_s){
		my ($s,$t) = (split /\t/,$po)[0,1];
		my $x1 = (1000*$s)/$conf_hash{"chr_length"};
		my $x2 = (1000*$t)/$conf_hash{"chr_length"};
		$svg->rect(x=>$x1+$forward_x,y=>$gene_y+$inter_gene_y,width=>$x2-$x1,height=>8,fill => $gcol{$genes_s{$po}},stroke=>$gcol{$genes_s{$po}},"stroke-width"=>0);
		$svg->text(x=>$x2+$forward_x,y=>$gene_y+$inter_gene_y,-cdata=>$genes_s{$po},"font-size"=>10,"font-family"=>"Arial","text-anchor"=>"right",fill=>"black");
		$inter_gene_y += 10;#每个基因的间隔是2，宽度是3
	}
	open OUT,">$out-graph.svg";
	print OUT $svg->xmlify;
	close OUT;
}

sub help_text {
	my $self = shift;
	return <<HELP
	USAGE:
	perl graph_main.pl line_svg -conf  -type  -output
HELP
}
1;

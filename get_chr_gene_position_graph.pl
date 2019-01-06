#! /usr/bin/perl

use utf8;
use strict;
use warnings;
use Getopt::Long;
use SVG;

my ($chr_len,$gene_posit_in,$out,$width,$heigth);

GetOptions(
	"c:s" => \$chr_len,
	"g:s" => \$gene_posit_in,
	"u:s" => \$out,
	"w:s" => \$width,
	"h:s" => \$heigth,
);

my $usage =<<"USAGE";
	Usage:  perl $0 -c chr_len -g gene_posit -u out [option]
		-c   chr_len
		-g   gene_posit
		-u   out
		-w   width
		-h   heigth
USAGE

if(!$chr_len or !$gene_posit_in){die $usage};

$width ||= 2000;
$heigth ||= 1000;
my $svg = SVG->new(width => $width,height => $heigth);

my @col = ("#EB756B","#F1A19A","#9FA01E","#B8BA4F","#28AD40","#1CB179","#6991CC","#94B3DF","#9A82BC","#B8A1CD","#D06BA5","#DD94BF","#00EE00","#00FFFF","#FFFF00","#32CD32","#8B658B");

my %genes_posit;
my %chr = &get_chr($chr_len,$gene_posit_in,\%genes_posit);

my $f_y = 100;
my $f_x = 100;
my $inter_x = 0;
&graph_chr(\%chr,$svg,$f_y,$f_x,$inter_x,\%genes_posit,\@col,$out);

sub get_chr{
	my ($in,$genes_in,$genes_posit) = @_;
	my %chr;
	open IN,$in or die $!;
	while(<IN>){
		chomp;
		my @t = split /\t/,$_;
		$chr{$t[0]} = $t[1];
	}
	close IN;
	open IN1,$genes_in or die $!;
	while(<IN1>){
		chomp;
		my @t = split /\t/,$_;
		push @{$genes_posit->{$t[0]}},[$t[1],$t[2],$t[3]];
	}
	close IN1;
	return %chr;
}

sub graph_chr{
	my ($chr,$svg,$f_y,$f_x,$inter_x,$genes_posit,$col,$out) = @_;
	my $n = 0;
	for my $i(0..9){
		my $n_y = $i*80;
		$svg->line(x1=>$f_x-20,y1=>$f_y+$n_y,x2=>$f_x-20,y2=>$f_y+$n_y+80,stroke=>"black",strokewidth=>0.5);
		$svg->line(x1=>$f_x-25,y1=>$f_y+$n_y,x2=>$f_x-20,y2=>$f_y+$n_y,stroke=>"black",strokewidth=>0.5);
		$svg->text(x=>$f_x-40,y=>$f_y+$n_y,-cdata=>"$i"."M","font-size"=>10,"font-family"=>"Arial","text-anchor"=>"right",fill=>"black");
	}
	for my $c (sort{$chr->{$b} <=> $chr->{$a}}keys %{$chr}){
		my $chr_l = ($chr->{$c}*800)/10000000;
		$svg->rect(x=>$f_x+$inter_x,y=>$f_y,width =>15,height =>$chr_l,fill =>$col->[$n],stroke=>$col->[$n],"stroke-width"=>1);
		for my $all (@{$genes_posit->{$c}}){
			my $genes_y = $f_y+($all->[0]*800)/10000000;
			$svg->line(x1=>$f_x+$inter_x+10,y1=>$genes_y,x2=>$f_x+$inter_x+20,y2=>$genes_y,stroke=>"black",strokewidth=>0.5);
			my $x = $f_x+$inter_x+15;
			my $y = $genes_y;
			$svg->text(x=>$f_x+$inter_x+25,y=>$genes_y+1,-cdata=>$all->[2],"font-size"=>10,"font-family"=>"Arial","text-anchor"=>"right",fill=>"black");
		}
		my $x = $f_x+$inter_x;
		my $y = $f_y-50;
		$svg->text(x=>$f_x+$inter_x,y=>$f_y-50,-cdata=>$c,transform=>"rotate(45,$x,$y)","font-size"=>10,"font-family"=>"Arial","text-anchor"=>"right",fill=>"black");
		$inter_x += 70;#(40 is intervect , 30 is width)
		$n++;
	}
	open OUT,">$out.svg";
	print OUT $svg->xmlify;
	close OUT;
}

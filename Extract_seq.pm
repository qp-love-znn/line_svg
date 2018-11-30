#!/usr/bin/perl -w

package Extract_seq;
require Exporter;
use strict;
use warnings;

our @ISA = qw(Exporter);
our @EXPORT = qw(set_file_tpye get_file_tpye get_gene_posit get_gene_seq get_gene_start_cds);
our @version = 1.0;

#sub get_gene_posit: get the gene position according to the gtf->my (%genes_posit) = $self->get_gene_posit($gtf)
#get_gene_seq : get the gene seq--->$self->get_gene_seq($fa,\%genes_posit)
#get_gene_start_cds:get the gene cds seq --->$self->get_gene_start_cds($gtf,$fa)

sub new{
	my $class = shift;
	my $self = {};
	my %parm = @_;
	bless $self;
	$self->{'file_tpye'} = $parm{'file_tpye'};
	return $self;
}

sub set_file_tpye{
	my $self = shift;
	my ($type) = @_;
	$self->{'file_tpye'} = $type;
	return $self->{'file_tpye'};
}

sub get_file_tpye{
	my $self = shift;
	if(!$self->{'file_tpye'}){
		print "file type have not begain\n";
	} else {
		print $self->{'file_tpye'}."\n";
	}
}
sub get_gene_posit{
	my $self = shift;
	my ($gtf) = @_;
#	if($type ne "gtf"){print "this not the gtf file\n";last};
	my %genes;
	open IN,$gtf or die $!;
	while(<IN>){
		chomp;
		my @t = split /\t/,$_;
		if($t[2] eq "exon"){
			my ($g) = $t[8] =~ /gene_id \"([^;]+)\";/;
			$genes{"strand"}{$g} = $t[6];
			if(exists $genes{$t[0]}{"start"}{$g}){
				if($genes{$t[0]}{"start"}{$g} > $t[3]){
					$genes{$t[0]}{"start"}{$g} = $t[3];
				}
				if($genes{$t[0]}{"end"}{$g} < $t[4]){
					$genes{$t[0]}{"end"}{$g} = $t[4];
				}
			} else {
				$genes{$t[0]}{"start"}{$g} = $t[3];
				$genes{$t[0]}{"end"}{$g} = $t[4];
			}
		}
	}
	close IN;
	return %genes;
}

sub get_gene_seq{
	my $self = shift;
	my ($fa,$start) = @_;
	open IN,$fa or die $!;
	$/=">";
	<IN>;
	while(<IN>){
		chomp;
		my @t = split /\n/,$_,2;
		my ($name) = $t[0] =~ /(^\S+)/;
		$t[1] =~ s/\n//g;
		if(exists $start->{$name}){
			for my $g (keys %{$start->{$name}->{"start"}}){
				my $gene_seq = substr($t[1],$start->{$name}->{"start"}->{$g}-1,$start->{$name}->{"end"}->{$g}-$start->{$name}->{"start"}->{$g});
				if($start->{"strand"}->{$g} eq "-"){
					$gene_seq =~ tr/ATGCatgc/TACGtacg/;
					print ">$g\n".reverse($gene_seq)."\n";
				} else{
					print ">$g\n$gene_seq\n";
				}
			}
		}
	}
	close IN;
}

sub get_gene_start_cds{
	my $self = shift;
	my ($gtf,$fa) = @_;
	open IN,$gtf or die $!;
	my (%genes_cds);
	while(<IN>){
		chomp;
		my @t = split /\t/,$_;
		my ($g) = $t[8] =~ /gene_id \"([^;]+)\";/;
		$genes_cds{"strand"}{$g} = $t[6];
		if(exists $genes_cds{$t[0]}{$g}){
			if($genes_cds{$t[0]}{$g} > $t[3]){
				$genes_cds{$t[0]}{$g} = $t[3];
			}
		} else {
			$genes_cds{$t[0]}{$g} = $t[3];
		}
	}
	close IN;
	open IN1,$fa or die $!;
	$/=">";
	<IN1>;
	while(<IN1>){
		chomp;
		my @t = split /\n/,$_,2;
		my ($name) = $t[0] =~ /(^\S+)/;
		$t[1] =~ s/\n//g;
		if(exists $genes_cds{$name}){
			for my $g (keys %{$genes_cds{$name}}){
				if($genes_cds{"strand"}{$g} eq "-"){
					my $st = $genes_cds{$name}{$g}-2000;
					my $seq = substr($t[1],$genes_cds{$name}{$g}-2001,2000);
					$seq =~ tr/ATGCatgc/TACGtacg/;
					print ">$g:$st-$genes_cds{$name}{$g}\n".reverse($seq)."\n";
				} else {
					my $st = $genes_cds{$name}{$g}-2000;
					my $seq = substr($t[1],$genes_cds{$name}{$g}-2001,2000);
					print ">$g:$st-$genes_cds{$name}{$g}\n$seq\n";
				}
			}
		}
	}
	close IN1;
}


1;

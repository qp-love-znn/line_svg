#!perl
use warnings;
use strict;

die "perl $0 <sorted bed>\n" if @ARGV != 1;


my %all;
my @arry;
my $chr_now;
my $stat_now;
my $end_now;
my $chr_old;
my $stat_old;
my $end_old;
open IN, $ARGV[0] or die $!;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	if($. ==1){
		$chr_old = $t[0];
		$stat_old = $t[1];
		$end_old = $t[2];
		push @arry,"$chr_old-$stat_old-$end_old";
	} else {
		$chr_now = $t[0];
		$stat_now = $t[1];
		$end_now = $t[2];
		if ($chr_old eq $chr_now){
			if($stat_now <= $end_old){
				if ($end_now >= $end_old){
					$end_old = $end_now;
					push @arry,"$chr_now-$stat_now-$end_now";
				} else {
					push @arry,"$chr_now-$stat_now-$end_now";
				}
			} else {
				if(@arry){
					print "$chr_old-$stat_old-$end_old\t".join("\t",@arry)."\n";
					@arry = ();
				} else {
					print "$chr_old-$stat_old-$end_old\n";
				}
				$chr_old = $chr_now;
				$stat_old = $stat_now;
				$end_old = $end_now;
			}
		} else {
			if(@arry){
				print "$chr_old-$stat_old-$end_old\t".join("\t",@arry)."\n";
				@arry = ();
			} else {
				print "$chr_old-$stat_old-$end_old\n";
			}
			$chr_old = $chr_now;
			$stat_old = $stat_now;
			$end_old = $end_now;
		}
	}
}
close IN;

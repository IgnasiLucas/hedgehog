#!/usr/bin/perl

 ##########################################################################
#										#	
#  Copyright (C) 2007-9 Jeffrey Ross-Ibarra <rossibarra@gmail.com>     	  	#
#                                                                         	#
#  This program is free software: you can redistribute it and/or modify   	#
#  it under the terms of the GNU General Public License as published by   	#
#  the Free Software Foundation, either version 3 of the License, or      	#
#  (at your option) any later version.                                    	#
#                                                                         	#
#  This program is distributed in the hope that it will be useful,        	#
#  but WITHOUT ANY WARRANTY; without even the implied warranty of         	#
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          	#
#  GNU General Public License <http://www.gnu.org/licenses/>			#
# 	for more details.   	                                             	#
#										#
 ##########################################################################

 ##########################################################################
#										#	
#  Thanks are due to R. Cartwright, P. Morrell, and K. Thornton    	        #
#										#
 ##########################################################################

use strict;
use Getopt::Std;

my $version="1.05";

# get options and input files
my %options=();
getopts("i:vcfmgntqsp",\%options) or help();

change() if $options{c};

# gets infiles
my %files = ();
$files{$options{i}} = 1 foreach(grep(-e, glob($options{i}))); 
my @files = keys(%files);
&help unless(@files);

foreach my $file (@files){

	open FILE, "<$file" or die help();

	my %alleles; # hash of polymorphic sites and their segregating bases 
	my @seqs=(); #AoA of all sequences and sites
	my @inds; # array of idnividual names from fasta file
	my %fourgam=(); # hash of recombo intervals $fourgam{$high}=$low
	my @polysites; # array of poylmorphic site numbers
	my @error; # array of sites with >2 segregating bases
	my %cuts=(); # hash of cutsites for nonrecombining bits
	my %gapped=(); #hash of sites with gaps
	my %missing=(); #hash of sites with missing data (N)
	my %trash=(); #hash of intervals to delete
	my %keep=(); #hash of intervals that gots recombination
	my $low=0; my $high=0; # positions for calculating rmin
	my $rmincount=0; # count of minimum recombo events
	my $yime=localtime();
	
	if ($options{v}){
		print "\n\nRMINCUTTER $version output for file $file on $yime\nOptions used: ";
		foreach my $used_opt ( keys(%options) ){ print " -", $used_opt unless $used_opt=~m/i/;  }
		print "\n\n";	
	}
	
	#turns sequences into AoA
	my $count=0; my $iterator=0; # iterators
	while(<FILE>){
		chomp $_;
		if($_=~m/^\>/){ 
			$inds[$count]=$_; 
			&format if $iterator != $count;
			$iterator=$count+1;
			next; 
		}
		last if $_=~/^$/;
		$_=~tr/[a-z]/[A-Z]/;
		$seqs[$count]=[split("",$_)];
		$count++;
	}
	close FILE;
	my @temp = @{ $seqs[0] }; 
	my $length=$#temp;
		
	#gets polymorphic sites
	for my $site (0.. $length){
		my %polycheck=();
		for my $i (0..$#seqs){ 	
		
			#skips gaps and missing data
			if($seqs[$i][$site]=~m/N/){ $missing{$site}=1; next }
			if($seqs[$i][$site]=~m/\-/){ $gapped{$site}=1; next }
			if($seqs[$i][$site]!~m/[A,C,T,G]/){ die "ERROR: $seqs[$i][$site] is not a recognized base\n\n"; }
			$polycheck{$seqs[$i][$site]}=1;
			
		}
		my @numsites=keys(%polycheck);
		$alleles{$site}=[@numsites];
		my $poly=$#numsites+1;
		if($poly>1){
			#skip missing (N)
			next if $options{n} && defined $missing{$site};
			next if $options{g} && defined $gapped{$site};
			next if $options{t} && $poly>2;
			push(@polysites,$site);
		}
	}
	
	@polysites=sort {$a <=> $b} @polysites;
	print "\nPAIRS OF SITES THAT FAIL 4-GAMETE TEST:\n\n" if $options{v};

	#does 4 gamete test
	#iterates over polymorphic sites only
	for my $psite (0..$#polysites){
		for my $nextsite ($psite+1..$#polysites){
			my %snps=();
			my $gametes=0; 
			
			# get all haplotype pairs in pop into hash
			for my $i (0..$#seqs){
			
				#ignore gaps and missing, but keep data from those columns
				if($seqs[$i][$polysites[$psite]] =~m/\-|N/ || $seqs[$i][$polysites[$nextsite]]=~m/\-|N/){ next; }
				
				$snps{ $seqs[$i][$polysites[$psite]] }->{ $seqs[$i][$polysites[$nextsite]] } = 1;
			}
			
			#get all alleles at both sites
			my @gam1=@{ $alleles{ $polysites[$psite] } }; 
			my @gam2=@{ $alleles{ $polysites[$nextsite] } }; 

			#count number of gametes for each pair of sites
			for my $a (0..$#gam1-1){
				for my $b ($a+1..$#gam1){
					for my $y (0..$#gam2-1){
						for my $z ($y+1..$#gam2){		
							$gametes=$snps{ $gam1[$a] }->{ $gam2[$y] } + $snps{ $gam1[$a] }->{ $gam2[$z] } + $snps{ $gam1[$b] }->{ $gam2[$y] } + $snps{ $gam1[$b] }->{ $gam2[$z] };						
						}
					}
				}
			}
	
			#make hash of intervals
			if($gametes>3){			
				print $polysites[$psite]+1, "\t", $polysites[$nextsite]+1, "\n" if $options{v};
				#if two intervals with same start, picks shortest
				if( defined $fourgam{ $polysites[$psite] } ){
					if( $polysites[$nextsite] < $fourgam{ $polysites[$psite] }){
						$fourgam{ $polysites[$psite] } = $polysites[$nextsite];
					}
				}
				else{ $fourgam{ $polysites[$psite] } = $polysites[$nextsite]; }
				#stop checking sites if find one that fails 4-gamete test
				last;
			}
		}
	}
	print "\n";
		
	my @starts=sort {$a<=>$b} keys(%fourgam);
	for my $start (0..$#starts){
		# -m gets minimum intervals ala Ross-Ibarra
		if($options{m}){
			if($start==0){
				$low=$starts[$start]; $high=$fourgam{$starts[$start]};			
			}
			elsif($starts[$start]>$low && $starts[$start]<$high){
				$low=$starts[$start];
				$high=$fourgam{ $starts[$start] } if $fourgam{$starts[$start]} < $high;
			}
			elsif($starts[$start]>=$high){
				$keep{ $low }=$high;
				$low=$starts[$start]; $high=$fourgam{ $starts[$start] };			
			}
			if($start==$#starts){ $keep{ $low } = $high; }
		}

		#throw out intervals following Hudson and Kaplan 1985	
		else{
			for my $twostart ( 0..$#starts ){		
				#don't self-compare or compare to already trashed intervals
				if( $twostart==$start || defined $trash{$starts[$twostart]} ){ next; }
				
				#step 1: ditch any segments that totally encompass other segments
				elsif($starts[$start]<=$starts[$twostart] && $fourgam{$starts[$start]}>=$fourgam{$starts[$twostart]}){
					$trash{$starts[$start]}=$fourgam{$starts[$start]};
					last;
				}
				#step 2: ditch any segments that overlap ( i1<i2<j1 in terms of Hudson and Kaplan 1985m )
				elsif($starts[$start]>=$starts[$twostart] && $starts[$start]<$fourgam{$starts[$twostart]}  ){
					$trash{$starts[$start]}=$fourgam{$starts[$start]};
					last;		
				}
			}
			#make hash of intervals to keep
			$keep{ $starts[$start] }=$fourgam{ $starts[$start] } unless defined $trash{ $starts[$start] };
		}
	}
	# determine cutpoints
	print "REGIONS WHERE RECOMBINATION IS ASSUMED TO HAVE OCCURED:\n\n" if $options{v};
	my $cutstart=0;
	foreach my $x ( sort {$a<=>$b} keys(%keep) ){
			print $x+1, "\t", $keep{$x}+1, "\n" if $options{v};
			if($x==$cutstart){
				$cutstart=$keep{$x};
				$rmincount++;
				next;
			}
			else{
				$cuts{ $cutstart } = $x;
				$cutstart=$keep{$x}; 
				$rmincount++;
			}
	}
	$cuts{ $cutstart }=$length if $cutstart < $length; 
	
	#rmin Hudson Kaplan
	print "\nRMIN: $rmincount\n" if $options{v};

	print "\nINDEPENDENT REGIONS TO BE WRITTEN TO OUTPUT:\n\n" if $options{v};
	my $piecemeal=1; # numbering for outfiles
	
	# print "greedy" regions that include all SNPSs ala Morrell
	if( $options{p} ){
		my $slicepoint=0;
		foreach my $x ( sort {$a<=>$b} keys(%keep) ){
			print "slice: $slicepoint\tx: $x\tkeepx: $keep{$x}\n";
			bprint($slicepoint, $keep{$x}-1, $piecemeal, $file, \@inds, \@seqs);
			$slicepoint=$keep{$x};
		}
		bprint($slicepoint, $length, $piecemeal, $file, \@inds, \@seqs);
	}
	
	# print normal intervals
	else{
		foreach my $x ( sort {$a<=>$b} keys(%cuts) ){
			bprint($x, $cuts{$x}, $piecemeal, $file, \@inds, \@seqs) 
		}
	}
	
	#report gapped (-) sites
	if($options{g}){
		print "\n\nSITES SKIPPED BECAUSE OF GAPS:\n\n";
		my @bob=sort {$a<=>$b} keys(%gapped);
		for my $site (@bob){
			print "$site\t"
		}
		print "none" if @bob==0;
		print "\n";
	}
	#report missing (N) sites
	if($options{n}){
		print "\n\nSITES SKIPPED BECAUSE OF MISSING DATA:\n\n";
		my @bob=sort {$a<=>$b} keys(%missing);
		for my $site (sort {$a<=>$b} keys(%missing)){
			print "$site\t"
		}
		print "none" if $#bob<1;
		print "\n";
	}
	#count and warn about sites with >2 bases
	$count=0; 
	foreach my $s (@polysites){ 
		my @temp=@{ $alleles{$s} }; 
		if( $#temp > 1){ 
			push(@error,$s);
			$count++;
		} 
	}
	if($count>0 && !$options{t} && !$options{s}){
		print STDERR "\nWARNING: In file $file the following sites have more than 2 bases";
		print STDERR "\nsegregating and should be checked as they may cause problems:\n\n";
		foreach my $err (sort {$a <=> $b} @error){
			print STDERR $err+1, " ( @{ $alleles{$err} } )\n";
		}
		print STDERR "\n";
	}
	#report sites skipped because >2 bases
	elsif($options{t}){
		print "\n\nSITES SKIPPED WITH >2 BASES SEGREGATING:\n\n";
		if($count>0){
			foreach my $err (sort {$a <=> $b} @error){
				print $err+1, " ( @{ $alleles{$err} } )\n";
			}
		}
		else{ print "none\n"; }
	}
	print "\n" if $options{v};
}

#SUBS##

sub help{
	print "\n********************************************\n*";
	print "\n* RminCutter $version Copyright 2007-9 Jeffrey Ross-Ibarra <rossibarra\@gmail.com>\n*";
	print "\n*\tRminCutter calculates rmin following the algorithm described in";
	print "\n*\tHudson, R. and N. Kaplan. 1985. Genetics 111:147-164\n*";	
	print "\n*\tRequired command line parameters:\n*\n*";
	print "\t\t-i 'filename' : input name of fasta file(s) (accepts * wildcard).\n*";
	print "\t\t\t\tOutput is printed to filename_1.out, filename_2.out, etc.\n*";
	print "\n*\tOptional command line parameters:\n*\n*";
	print "\t\t-v : for verbose output\n*";
	print "\t\t-q : suppress top header line on output files\n*";
	print "\t\t-s : suppress warning to STDERR\n*";
	print "\t\t-f : add '.fasta' to the end of output files for easier manipulation by\n*"; 
	print "\t\t     other programs\n*";
	print "\t\t-m : for smaller set of recombination intervals\n*";
	print "\t\t     This option forces reduction of recombination\n*";
	print "\t\t     intervals to the smallest range of sequence\n*";
	print "\t\t     which can explain the observed data.\n*";
	print "\t\t     It does NOT give regions in which recombination DID occur.\n*";
	print "\t\t-g : skip sites with gaps\n*";
	print "\t\t-p : uses a 'greedy' algorithm which uses all SNPs\n*";
	print "\t\t     As a result, each region is likely to still encompass a recombination\n*";
	print "\t\t     event and thus violate the assumptions of many models.\n*";
	print "\t\t-n : skip sites with missing data\n*";
	print "\t\t-t : skip sites with >2 bases segregating data\n*";
	print "\t\t-c : changelog since version 1.0\n*";
	print "\n********************************************\n";
	die "\n";
}

sub change{
	print "\n********************************************\n*";
	print "\n* RminCutter $version Copyright 2009 Jeffrey Ross-Ibarra <rossibarra\@gmail.com>\n*";
	print "\n*\tChanges since v 1.0:\n*";
	print "\n*\t12/03/07: blank lines at end of fasta file no longer cause bug";
	print "\n*\t12/04/07: added ability to ignore sites with >2 bases segregating";
	print "\n*\t12/04/07: fixed bug in missing/gap order of removal";
	print "\n*\t12/04/07: checked w/ K. Thornton's compute: using -t -g gives identical results";
	print "\n*\t02/23/08: fixed verbose output of gapped sites";
	print "\n*\t02/23/08: added the 'greedy' algorithm ala Morrell";
	print "\n*\t02/24/08: added format and symbol checking";
	print "\n*\t02/24/08: added silent option to skip warnings";
	print "\n*\t02/24/08: added bprint subroutine";
	print "\n*\t02/24/08: changed verbose output: text for minimum regions, options used";
	print "\n*\t12/27/08: fixed bug in 'greedy' alogrithm that wrote region end SNPs twice";
	print "\n*\t06/28/09: fixed bug in 'greedy' alogrithm that made some printed regions fail 4-gamete test";
	print "\n*\t06/28/09: minor changes to file opening\n*\n*";
	print "\n********************************************\n";
	die "\n";
}

sub format{
	print "\nUnexpected formatting found.";
	print "\nPlease format your fasta file as follows, with each sequence all on one line:";
	print "\n\n>name1\nTTCGATCGATCGATCG\n>name2\nATC-ATTGATCGATCG\n>name3\nATCGATNGATCGATCC\n";
	die "\n";
}

# takes start, end, iterator, filename, inds array, seqs array
sub bprint{
	my $start=$_[0];
	my $end=$_[1];
	my $file=$_[3];
	my @inds = @{$_[4]};
	my @seqs = @{$_[5]};
	
	print $start+1, "\t", $end+1, "\n" if $options{v};

	my $outfile;
 	if ($options{f}){ $outfile = $file . "_" . $_[2] . ".fasta"; }
	else{ $outfile = $file . "_" . $_[2] . ".out"; }
	open OUT, ">$outfile";
	unless( $options{q} ) { 
		print OUT "# RminCutter $version output from $file regions ", $start+1, " to  ", $end+1, "\n";
	}
	for my $i (0..$#seqs){
		print OUT "$inds[$i]\n";
		for my $bp ($start..$end){
			print OUT $seqs[$i]->[$bp];
		}
		print OUT "\n"
	}
	close OUT;
	$_[2]++;	
}

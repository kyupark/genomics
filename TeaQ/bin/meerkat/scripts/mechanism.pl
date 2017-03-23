#!/usr/bin/perl
# mechanism.pl
# assign mechanism to each variant
# Author: Lixing Yang, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
# Email: lixing_yang@hms.harvard.edu
# Input: sorted bam file, repeat file
# Output:
# TEI: transposable element insertion, single or multiple TE insertion.
# TEA: Alternative TE, usually a deletion with insertion event is called, deletion is a TE in reference genome and insertion is the same type of TE in reference genome. It should be arisen from sequence divergence of TE.
# VNTR: variable number of tandem repeat, deletion or insertion of satellite repeat, simple repeat or low complexity repeat.
# NAHR: non-allellic homologous recombination, >100bp homology.
# alt-EJ: alternative end joining, 3-100bp homology.
# NHEJ: non-homologous end joining, 0-2bp homology or 1-10bp insertion at deletion break point.
# FoSTeS: fork stalling and template switching, template switch, >10bp insertion at deletion break points.
# NA: unclassified.
# i.e. perl scripts/mechanism.pl -R /db/hg18/rmsk-hg18.txt -b

use strict;
use Getopt::Std;
use FindBin '$Bin';
my $version = 'v.0.174';

my %opts = (D=>"null", o=>1, T=>100000, z=>1000000000);
getopts("D:b:o:z:R:T:t:z:h", \%opts);
my $tmp_dir = $opts{D};
my $bamfile = $opts{b};
my $include_other = $opts{o};
my $rmskfile = $opts{R};
my $te_size_max = $opts{T};
my $threads_bwa = $opts{t};
my $sv_size_cutoff = $opts{z};
my $del_ins_size_cutoff_d = 0.8; # size ratio of del and ins in del_ins events
my $del_ins_size_cutoff_u = 1.2;
my $ovl = 0.8; # overlap of a predicted events and a annotated TE

&print_usage unless (defined($opts{b}) and defined($opts{R}));
&print_usage if (defined($opts{h}));

$Bin =~ /scripts/;
my $bin_path          = $` . 'bin/';
my $bamreader_command = $bin_path . 'bamreader';

my $time0 = time;
my $local0 = localtime($time0);
print STDERR "$local0 Mechanism $version started\n";

my $prefix = $bamfile;
if ($prefix =~ /.bam$/)
{
	$prefix = $`;
}
if ($prefix =~ /.sorted$/)
{
	$prefix = $`;
}
(my $working_prefix = $prefix) =~ s{.*/}{};
my $is_actual_dir = $prefix;
my $is_files_dir;
if("null" ne $tmp_dir) {
	if(-e "$tmp_dir/$working_prefix") {
		$is_actual_dir = "$tmp_dir/$working_prefix";
		$is_files_dir = $is_actual_dir."/".$working_prefix;
	}
}
my $intra_refine_type = $prefix.'.intra.refined.typ.sorted';
my $inter_refine_type = $prefix.'.inter.refined.typ.sorted';
my $outfile = $prefix.'.variants';

if(-e $is_files_dir) {
$intra_refine_type = $is_files_dir.'.intra.refined.typ.sorted';
$inter_refine_type = $is_files_dir.'.inter.refined.typ.sorted';
$outfile = $is_files_dir.'.variants';
	
}
if("null" ne $tmp_dir) {
	print("$bamreader_command mechanism -t $threads_bwa -x $prefix -D $is_actual_dir -o $include_other -R $rmskfile -T $te_size_max -sv $sv_size_cutoff\n");
	system("$bamreader_command mechanism -t $threads_bwa -x $prefix -D $is_actual_dir -o $include_other -R $rmskfile -T $te_size_max -sv $sv_size_cutoff");
} else {
	print("$bamreader_command mechanism -t $threads_bwa -x $prefix -o $include_other -R $rmskfile -T $te_size_max -sv $sv_size_cutoff\n");
	system("$bamreader_command mechanism -t $threads_bwa -x $prefix -o $include_other -R $rmskfile -T $te_size_max -sv $sv_size_cutoff");
}

#use Cwd qw();
#my $cur_path = Cwd::abs_path();
#my $target_path = "$cur_path/called/$working_prefix"; 
#use File::Path qw/make_path/;
#make_path($target_path);
#system("mv $is_actual_dir/* $target_path");

my $time1 = time;
my $local1 = localtime($time1);
my $difference = $time1 - $time0;
my $seconds    =  $difference % 60;
$difference = ($difference - $seconds) / 60;
my $minutes    =  $difference % 60;
$difference = ($difference - $minutes) / 60;
print STDERR "$local1 Finished\n";
print STDERR "Time used: $difference:$minutes:$seconds\n";


sub print_usage
{
	die "mechanism.pl [options]
	-b FILE	sorted and indexed bam file, required
	-o INT	[0/1], include rmsk type \"Other\" in TE, default 1
	-t INT	max size of TE, default 100,000
	-z INT	size limit of SVs to be processed, default 1,000,000,000
	-R STR	/path/to/repeat_mask_file, required, can be downloaded from UCSC
	-h help\n";
}

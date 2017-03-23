#!/usr/bin/perl
# pre_process.pl
# generate blacklist of positions, insert size distribution and unmapped reads
# Author: Lixing Yang, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
# Email: lixing_yang@hms.harvard.edu
# Input: sorted bam file, coverage cutoff, insertsize range for plotting
# i.e. perl scripts/pre_process.pl -I /db/hg18/hg18_bwa_idx/hg18.fasta -A /db/hg18/hg18.fasta.fai -W /opt/bwa/ -b
# Output: xxx.blacklist (positions to be ignored), xxx.unmapped.fq.gz (unmapped reads in fastq format), xxx.isinfo (insert size distribution), xxx.pdf (insert size distribution plot)
# System requirement:
# samtools 0.1.5 or above
# R 2.6.1 or above

use strict;
use Getopt::Std;
use FindBin '$Bin';
my $version = 'v.0.189';

my %opts = (k=>500, r=>1000, n=>10000, l=>1, t=>1, q=>15, c=>35, s=>20, u=>0, f=>100, N=>5, R=>'null', P=>'all');
getopts("A:b:c:C:f:I:k:l:Mn:N:P:q:r:R:s:S:t:u:W:h", \%opts);
my $tmp_dir = "tmp";
my $reference_fai = $opts{A};
my $bamfile = $opts{b};
my $clip_cut = $opts{c};
my $caldepth_path = $opts{C};

my $n_alter_map = $opts{f};
my $bwa_index = $opts{I};
my $blacklist_cutoff = $opts{k};
my $remap_cl = $opts{l};
my $is_creating_um_sc = 1 if(defined($opts{M}));
my $number_reads = $opts{n};
my $numberofn = $opts{N};

my $step = $opts{P};
my $qual_trim = $opts{q};
my $plotrange = $opts{r};
my $black_rg = $opts{R};
my $sr_cut = $opts{s};
my $samtools_path = $opts{S};
my $threads_bwa = $opts{t};
my $uupair = $opts{u};
my $bwa_path = $opts{W};

my $samtools_command;
if (defined($opts{S}))
{
	$samtools_path .= '/' unless ($samtools_path =~ /\/$/);
	$samtools_command = $samtools_path.'samtools';
}
else
{
	$samtools_command = 'samtools';
}
my $bwa_command;
if (defined($opts{W}))
{
	$bwa_path .= '/' unless ($bwa_path =~ /\/$/);
	$bwa_command = $bwa_path.'bwa';
}
else
{
	$bwa_command = 'bwa';
}
$Bin =~ /scripts/;
my $bin_path = $`.'bin/';
my $bamreader_command = $bin_path.'bamreader';

&print_usage unless (defined($bamfile));
&print_usage if (defined($opts{h}));
my $time0 = time;
my $local0 = localtime($time0);
print STDERR "$local0 Pre-process $version started\n";

my $is = 1 if ($step eq 'all' or $step eq 'is');
my $cl = 1 if ($step eq 'all' or $step =~ /cl/);
my $cls =        0;
$cls =           1 if ($step eq 'cl1');
$cls =           2 if ($step eq 'cl2');

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
use File::Path qw/make_path/;
$is_actual_dir = "$tmp_dir/$working_prefix";
$is_files_dir = $is_actual_dir."/".$working_prefix;
make_path("$is_files_dir");
#print("$is_files_dir\n");

my $unmappedfile = $prefix.'.unmapped.fq.gz';
my $scfile = $prefix.'.softclips.fq.gz';
if(-e $is_actual_dir) {
	$unmappedfile = $is_actual_dir.'.unmapped.fq.gz';
	$scfile = $is_actual_dir.'.softclips.fq.gz';
}

my %blackrg;
if (-e $black_rg)
{
	open FILE, "<$black_rg";
	my $newline;
	while ($newline = <FILE>)
	{
		chomp $newline;
		$blackrg{$newline} = 1;
	}
	close FILE;
}

my $rdist_unmapfile = $prefix.'.unmapped.rdist';
my $rdist_scfile = $prefix.'.softclips.rdist';
my $rfile = $prefix.'.insert.r';
my $pdffile = $prefix.'.pdf';
my $cldir = $prefix.'/';
my $resultfile = $prefix.'.isinfo';
my $logfile = $prefix.'.pre.log';
if(-e $is_actual_dir) {
	$rdist_unmapfile = $is_files_dir.'.unmapped.rdist';
	$rdist_scfile = $is_files_dir.'.softclips.rdist';
	$rfile = $is_files_dir.'.insert.r';
	$pdffile = $is_files_dir.'.pdf';
	$cldir = $is_files_dir.'/';
	$resultfile = $is_files_dir.'.isinfo';
	$logfile = $is_files_dir.'.pre.log';
}
	
if ($is) {
	my $preprocess_cmd = "$bamreader_command preprocess -t $threads_bwa -D $is_actual_dir -F $bwa_index -l $remap_cl -b $blacklist_cutoff -x $prefix -r $black_rg -c $clip_cut -s $sr_cut -n $number_reads -q $qual_trim -N $numberofn -f $n_alter_map ";
	if ($uupair) {
		$preprocess_cmd = $preprocess_cmd."-u ";
		if($is_creating_um_sc) {
			$preprocess_cmd = $preprocess_cmd."-M ";
		}		
		$preprocess_cmd = $preprocess_cmd."$bamfile > $logfile"; 
		print("$preprocess_cmd\n");
		system("$preprocess_cmd");
	} else {
		if($is_creating_um_sc) {
			$preprocess_cmd = $preprocess_cmd."-M ";
		}
		$preprocess_cmd = $preprocess_cmd."$bamfile > $logfile";
		print("$preprocess_cmd\n");
		system("$preprocess_cmd");
	}
#my @rg;

#if(-e $is_actual_dir) {
#	opendir(DIR, $is_files_dir) || die "Can't open $is_files_dir";
#	# Read the dir ignoring . and ..
#	my @files = grep !/^\.\.?$/, readdir DIR ;
#	foreach my $infile (@files)
#	{
#		if ($infile =~ /.is$/)
#		{
#			my $filesize = -s $is_files_dir.'/'.$infile;
#			next unless ($filesize);
#			push @rg, $`;
#		}
#	}
#} else {
#	opendir(DIR, $prefix) || die "Can't open $prefix";
#	# Read the dir ignoring . and ..
#	my @files = grep !/^\.\.?$/, readdir DIR ;
#	foreach my $infile (@files)
#	{
#		if ($infile =~ /.is$/)
#		{
#			my $filesize = -s $prefix.'/'.$infile;
#			next unless ($filesize);
#			push @rg, $`;
#		}
#	}
#}

#open RFILE, ">$rfile";
#print RFILE "pdf(\"$pdffile\")\n";
#my $rgi=0;
#foreach (@rg)
#{
#	next if ($blackrg{$_});
#	my $insertsizefile = $prefix.'/'.$_.'.is';
#	if(-e $is_actual_dir) {
#		$insertsizefile = $is_files_dir.'/'.$_.'.is';
#	}
#	my $rg = $_;
#	my $rgname = 'rg'.$rgi;
#	my $filesize = -s $insertsizefile;
#	if ($filesize > 0)
#	{
#		print RFILE "$rgname=read.table(\"$insertsizefile\")\n";
#		if ($rgi)
#		{
#			print RFILE "data=c(data,$rgname\$V1)\n";
#		}
#		else
#		{
#			print RFILE "data=$rgname\$V1\n";
#		}
#		print RFILE "m=mean($rgname\[,1\])
#n=median($rgname\[,1\])
#sd=sd($rgname\[,1\])
#write(\"Mean insert size:\t$rg\", append=T, file=\"$resultfile\")
#write(m, append=T, file=\"$resultfile\")
#write(\"Median insert size:\t$rg\", append=T, file=\"$resultfile\")
#write(n, append=T, file=\"$resultfile\")
#write(\"Standard deviation of insert size:\t$rg\", append=T, file=\"$resultfile\")
#write(sd, append=T, file=\"$resultfile\")\n";
#	}
#	$rgi++;
#}
#print RFILE "hist(data,xlim=c(0,$plotrange),breaks=1000,xlab=\"Insert size\",main=\"\")
#unmap=read.table(\"$rdist_unmapfile\")
#sc=read.table(\"$rdist_scfile\")
#barplot(unmap[,3],col=1,space=0,names.arg=unmap[,1],main=\"Read length distribution of unmapped reads\",xlab=\"Read length (bp)\",ylab=\"Fraction\")
#barplot(sc[,3],col=1,space=0,names.arg=sc[,1],main=\"Read length distribution of soft clipped reads\",xlab=\"Read length (bp)\",ylab=\"Fraction\")
#";
#close RFILE;
#print("Rscript $rfile");
#system("Rscript $rfile");
#system "rm $rfile";
#if(-e $is_actual_dir) {
#foreach (@rg)
#{
#	next if ($blackrg{$_});
#	my $insertsizefile = $is_actual_dir.'.is.'.$_;
##	system "rm $insertsizefile";
#}
#} else {
#foreach (@rg)
#{
#	next if ($blackrg{$_});
#	my $insertsizefile = $prefix.'.is.'.$_;
#	if(-e $is_actual_dir) {
#		$insertsizefile = $is_actual_dir.'.is.'.$_;
#	}
##	system "rm $insertsizefile";
#}
#}
}

#my @rg;
#if(-e $is_actual_dir) {
#	opendir(DIR, $is_files_dir) || die "Can't open $prefix";
## Read the dir ignoring . and ..
#	my @files = grep !/^\.\.?$/, readdir DIR ;
#	foreach my $infile (@files)
#{
#	if ($infile =~ /.is$/)
#	{
##		print("$infile\n");
#		my $filesize = -s $is_files_dir.'/'.$infile;
#		next unless ($filesize);
#		my $rgname = $`;
##		print("$rgname\n");
#		my $cl_fq1 = $cldir.$rgname.'_1.fq';
##		print("$cl_fq1\n");
#		next unless (-s $cl_fq1 > 0);
#		push @rg, $rgname;
#	}
#}
#} else {
#	opendir(DIR, $prefix) || die "Can't open $prefix";
## Read the dir ignoring . and ..
#	my @files = grep !/^\.\.?$/, readdir DIR ;
#	foreach my $infile (@files)
#{
#	if ($infile =~ /.is$/)
#	{
##		print("$infile\n");
#		my $filesize = -s $prefix.'/'.$infile;
#		if(-e $is_actual_dir) {
#			$filesize = -s $is_actual_dir.'/'.$infile;
#		}
#		next unless ($filesize);
#		my $rgname = $`;
#		my $cl_fq1 = $cldir.$rgname.'_1.fq';
#		next unless (-s $cl_fq1 > 0);
#		push @rg, $rgname;
#	}
#}
#}
#if ($cl or $remap_cl)
#{
#if ($rg[0] eq 'none')
#{
#	my $cl_fq1 = $cldir.'none_1.fq';
#	my $cl_fq2 = $cldir.'none_2.fq';
#	my $cl_sai1 = $prefix.'.cl_1.sai';
#	my $cl_sai2 = $prefix.'.cl_2.sai';
#	my $cl_bam = $prefix.'.cl.bam';
#	my $cl_sort = $prefix.'.cl.sorted';
#	my $cl_sortbam = $prefix.'.cl.sorted.bam';
#	if(-e $is_actual_dir) {
#		$cl_sai1 = $is_files_dir.'.cl_1.sai';
#		$cl_sai2 = $is_files_dir.'.cl_2.sai';
#		$cl_bam = $is_files_dir.'.cl.bam';
#		$cl_sort = $is_files_dir.'.cl.sorted';
#		$cl_sortbam = $is_files_dir.'.cl.sorted.bam';
#	}
#	if (-s $cl_fq1 > 0)
#	{
#		goto CLS1 if ($cls == 2);
#		system "$bwa_command aln $bwa_index -l 40 -k 2 -t $threads_bwa $cl_fq1 > $cl_sai1  2>>bwa.err";
#		system "$bwa_command aln $bwa_index -l 40 -k 2 -t $threads_bwa $cl_fq2 > $cl_sai2  2>>bwa.err";
#		goto DONE if ($cls == 1);
#CLS1:
#		print("$bwa_command sampe -P -N $n_alter_map $bwa_index $cl_sai1 $cl_sai2 $cl_fq1 $cl_fq2 2>>bwa.err | $samtools_command view -bt $reference_fai -o $cl_bam -\n");	
#		system("$bwa_command sampe -P -N $n_alter_map $bwa_index $cl_sai1 $cl_sai2 $cl_fq1 $cl_fq2 2>>bwa.err | $samtools_command view -bt $reference_fai -o $cl_bam -");
#		print("sambamba sort -t $threads_bwa -o $cl_sortbam $cl_bam\n");
#		system("sambamba sort -t $threads_bwa -o $cl_sortbam $cl_bam");
#		print("$samtools_command index $cl_sortbam\n");
#		system("$samtools_command index $cl_sortbam");
#		#system "rm $cl_sai1";
#		#system "rm $cl_sai2";
##		system "rm $cl_bam";
#	}
#}
#else
#{
#	my @bam;
#	foreach (@rg)
#	{
#		next if ($blackrg{$_});
#		my $rgname = $_;
#		$rgname =~ s/\'/\\\'/g;
#		my $cl_fq1 = $cldir.$rgname.'_1.fq';
#		my $cl_fq2 = $cldir.$rgname.'_2.fq';
#		my $cl_sai1 = $cldir.$rgname.'_1.sai';
#		my $cl_sai2 = $cldir.$rgname.'_2.sai';
#		my $cl_bam = $cldir.$rgname.'.bam';
#		push @bam, $cl_bam;
#		goto CLS2 if ($cls == 2);
#		print("$bwa_command aln $bwa_index -l 40 -k 2 -t $threads_bwa $cl_fq1 > $cl_sai1  2>>bwa.err\n");
#		system("$bwa_command aln $bwa_index -l 40 -k 2 -t $threads_bwa $cl_fq1 > $cl_sai1  2>>bwa.err");
#		print("$bwa_command aln $bwa_index -l 40 -k 2 -t $threads_bwa $cl_fq2 > $cl_sai2  2>>bwa.err\n");
#		system("$bwa_command aln $bwa_index -l 40 -k 2 -t $threads_bwa $cl_fq2 > $cl_sai2  2>>bwa.err");
#		next if ($cls == 1);
#CLS2:	
#print("$bwa_command sampe -P -N $n_alter_map $bwa_index $cl_sai1 $cl_sai2 $cl_fq1 $cl_fq2 2>>bwa.err |perl -e 'while (<>){chomp;\@a=split(/\\t/,\$_);if(\$a[1] =~ /:/){print \"\$_\\n\";}else{if (\$a[11]){\$a[11]=\"RG:Z:$_\\t\".\$a[11];}else{\$a[11]=\"RG:Z:$_\";}print join(\"\\t\",\@a),\"\\n\";}}' | $samtools_command view -bt $reference_fai -o $cl_bam -\n");
#system("$bwa_command sampe -P -N $n_alter_map $bwa_index $cl_sai1 $cl_sai2 $cl_fq1 $cl_fq2 2>>bwa.err |perl -e 'while (<>){chomp;\@a=split(/\\t/,\$_);if(\$a[1] =~ /:/){print \"\$_\\n\";}else{if (\$a[11]){\$a[11]=\"RG:Z:$_\\t\".\$a[11];}else{\$a[11]=\"RG:Z:$_\";}print join(\"\\t\",\@a),\"\\n\";}}' | $samtools_command view -bt $reference_fai -o $cl_bam -");
#		#system "rm $cl_sai1";
#		#system "rm $cl_sai2";
#	}
#	goto DONE if ($cls == 1);
#	my $bam_merge = $prefix.'.cl.bam';
#	my $cl_sort = $prefix.'.cl.sorted';
#	my $cl_sortbam = $prefix.'.cl.sorted.bam';
#	if(-e $is_actual_dir) {
#		$bam_merge = $is_files_dir.'.cl.bam';
#		$cl_sort = $is_files_dir.'.cl.sorted';
#		$cl_sortbam = $is_files_dir.'.cl.sorted.bam';
#	}
#	if (@bam > 1)
#	{
#		print("sambamba merge -t $threads_bwa -l 1 $bam_merge @bam\n");
#		system("sambamba merge -t $threads_bwa -l 1 $bam_merge @bam");
#	}
#	else
#	{
#		print("cp $bam[0] $bam_merge\n");
#		system("cp $bam[0] $bam_merge");
#	}
#	print("sambamba sort -t $threads_bwa -o $cl_sortbam $bam_merge\n");
#	system("sambamba sort -t $threads_bwa -o $cl_sortbam $bam_merge");
#	print("$samtools_command index $cl_sortbam\n");
#	system("$samtools_command index $cl_sortbam");
#	foreach (@rg)
#	{
#		next if ($blackrg{$_});
#		my $cl_bam = $cldir.$_.'.bam';
#		my $rgname = $_;
#		$rgname =~ s/\'/\\\'/g;
#		if (-e $cl_bam)
#		{
#			my $cl_bam = $cldir.$rgname.'.bam';
##			system "rm $cl_bam";
#		}
#	}
##	system "rm $bam_merge";
#}
#}

DONE: my $time1 = time;
my $local1 = localtime($time1);
my $difference = $time1 - $time0;
my $seconds    =  $difference % 60;
$difference = ($difference - $seconds) / 60;
my $minutes    =  $difference % 60;
$difference = ($difference - $minutes) / 60;
print STDERR "$local1 Pre-process finished\n";
print STDERR "Time used: $difference:$minutes:$seconds\n";

sub print_usage
{
	die "pre_process.pl [options]
	-b FILE	sorted and indexed bam file, required
	-k INT	min coverage required for a nucleotide to be ignored (ignore position that covered by too many reads), 0 turn this function off, default 500
	-r INT	range of insert size to be plotted, insert size larger than such won't be used calculating distribution, default 1000
	-n INT	number of reads per read group to be used calculating insert size distribution, 0 use all reads, default 10000
	-l INT	[0/1], re-map read pairs of softclipped reads and unmapped reads, default 1
	-q INT	Read trimming parameter. Equivalent to BWA's -q option, default 15
	-c INT	bp to be cut off from beginning and end of unmapped and soft clipped reads, such reads were used to find smaller events, default 35
	-s INT	bp to be cut off from beginning and end of unmapped and soft clipped reads, such reads were used in split reads mapping, must be same as -s in meerkat.pl, default 20
	-u INT	[0/1], process uu pair (both reads unmapped in a pair) into 4 pairs, default 0
	-f INT	number of alternative mappings to print in XA tag for clipped alignments, default 100
	-N INT	Clipped reads and split reads must have <= INT Ns, default 5
	-t INT	number of threads used in bwa alignment, default 1
	-R FILE	file name of read group to be ignored, one read group ID per line
	-I STR	/path/to/reference/bwa_index for bwa alignment, required if -l 1
	-A STR	/path/to/reference/fasta.fai for bwa alignment, required if -l 1
	-S STR	/path/to/samtools, path only, not the command, no need to specify if samtools is in PATH
	-W STR	/path/to/bwa, path only, not the command, no need to specify if bwa is in PATH
	-P STR	specify step to run, [all|is|cl], default all
		is: extract unmapped, soft clipped reads, calculate insert size distribution
		cl: map soft clipped read pairs to reference genome
		all: run all above steps
	-h help\n";
}

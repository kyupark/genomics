#!/usr/bin/perl
# meerkat.pl
# Identify SVs by read pairs and split reads, give precise break points by local alignment
# Author: Lixing Yang, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
# Email: lixing_yang@hms.harvard.edu
# System requirements:
# samtools 0.1.5 or above
# BWA 0.5.7 or above
# NCBI blast 2.2.10 or above
# i.e. perl scripts/meerkat.pl -F /db/hg18/hg18_fasta/ -W /opt/bwa/ -B /opt/blast/bin/ -b
# Revised by Euncheon Lim, since Jun 17, 2016, euncheon_lim@hms.harvard.edu
use strict;
no warnings;
use Getopt::Std;
#use Bio::DB::Fasta;
use FindBin '$Bin';
my $version = 'v.0.189';
my %opts    = (
	k => 1,
	d => 3,
	D => "null",
	p => 2,
	o => 0,
	q => 1,
	z => 1000000000,
	s => 20,
	m => 1,
	a => 1,
	u => 0,
	Q => 0,
	l => 1,
	t => 1,
	R => 'null',
	P => 'all'
);
getopts("b:k:d:c:D:p:o:q:z:s:m:a:u:Q:g:f:l:t:R:F:S:W:B:P:n:h", \%opts);
my $ad_align         = $opts{a};
my $bamfile          = $opts{b};
my $blastall_path    = $opts{B};
my $sd_cutoff_cl     = $opts{c};
my $sd_cutoff_disc   = $opts{d};
my $tmp_dir			 = $opts{D};
my $alt_map_max_clip = $opts{f};
my $reference_path   = $opts{F};
my $alt_map_max      = $opts{g};
my $blacklist        = $opts{k};
my $clip             = $opts{l};
my $remove_dup       = $opts{m};
my $line             = $opts{n};
my $support_mpf      = $opts{o};
my $support_mps      = $opts{p};
my $step             = $opts{P};
my $support_reads	 = $opts{q};
my $min_mapq         = $opts{Q};
my $black_rg         = $opts{R};
my $cut_sr           = $opts{s};
my $samtools_path    = $opts{S};
my $threads_bwa      = $opts{t};
my $use_all_align    = $opts{u};
my $bwa_path         = $opts{W};
my $sv_size_cutoff   = $opts{z};

$alt_map_max = -1 unless (defined($opts{g}));
$alt_map_max_clip = -1 unless (defined($opts{f}));
$ad_align = 0 if $use_all_align;
$sd_cutoff_cl = $sd_cutoff_disc unless (defined($opts{c}));

if (defined($opts{h})) {
	&print_usage;
	die;
}
die "bam file not specified\n" unless (defined($bamfile));
my $prefix = $bamfile;
if ($prefix =~ /.bam$/) {
	$prefix = $`;
}
#if ($prefix =~ /.sorted$/) {
#	$prefix = $`;
#}


#	control progress
my $discordant = 1
  if ($step eq 'all' or $step eq 'dc');    # extract discordant read pairs
my $cluster = 1
  if ($step eq 'all' or $step =~ /cl/);    # generate discordant read pair clusters
my $mpd = 1
  if ($step eq 'all' or $step eq 'mpd');    # call candidates from discordant clusters
my $alg = 1
  if ($step eq 'all' or $step =~ /alg/);    # construct search space based on mpd candidates and align split reads to search space
my $srd = 1
  if ($step eq 'all' or $step eq 'srd');    # call SVs from pseudo split reads
my $filter = 1
  if ($step eq 'all' or $step eq 'srd' or $step eq 'ft');    # filter outputs
my $blast_bp = 1
  if ($step eq 'all' or $step eq 'rf');                      # blast break points reads to break points regions
my $cl = 0;
$cl = 1 if ($step eq 'cl1');
$cl = 2 if ($step eq 'cl2');
$cl = 3 if ($step eq 'cl3');
my $algs = 0;
$algs = 1 if ($step eq 'alg1');
$algs = 2 if ($step eq 'alg2');
die "path to reference fasta files not specified\n"
  if ($alg and !(defined($opts{F})));
my $prefix = substr($bamfile, 0, -4);

#if ($prefix =~ /.sorted$/) {
#	$prefix = $`;
#}

(my $working_prefix = $prefix) =~ s{.*/}{};
my $is_actual_dir = $prefix;
my $is_files_dir;
if("null" ne $tmp_dir) {
	if(-e "$tmp_dir/$working_prefix") {
		$is_actual_dir = "$tmp_dir/$working_prefix";
		$is_files_dir = $is_actual_dir."/".$working_prefix;
	}
}
if (defined($step)) {
	die "wrong step specified\nmust be one of the following: all|dc|cl|mpd|alg|srd|rf, default all\n"
	unless ($step =~ /all|dc|cl|mpd|alg|srd|ft|rf/);
}
my $blacklistfile = $prefix . '.blacklist.gz' if ($opts{k});
if(-e $is_files_dir) {
	$blacklistfile = $is_files_dir.'.blacklist.gz' if ($opts{k});
}
my ($newline, %is);

# parse median and standard deviation of insert size from file
my $isinfofile = $prefix . '.isinfo';
if(-e $is_files_dir) {
	$isinfofile = $is_files_dir. '.isinfo';
}
die "$isinfofile file not exist\n" unless (-e $isinfofile);
open FILE, "<$isinfofile";
while ($newline = <FILE>) {
	chomp $newline;
	if ($newline =~ /Read length/) {
		my ($trash, $rg) = split("\t", $newline);
		$newline = <FILE>;
		chomp $newline;
		$is{$rg}{'rl'} = $newline;
		if ($is{'rlu'}) {
			$is{'rlu'} = $is{$rg}{'rl'} if ($is{$rg}{'rl'} > $is{'rlu'});
		} else {
			$is{'rlu'} = $is{$rg}{'rl'};
		}
		if ($is{'rld'}) {
			$is{'rld'} = $is{$rg}{'rl'} if ($is{$rg}{'rl'} < $is{'rld'});
		} else {
			$is{'rld'} = $is{$rg}{'rl'};
		}
	}
	if ($newline =~ /Median/) {
		my ($trash, $rg) = split("\t", $newline);
		$newline = <FILE>;
		chomp $newline;
		$is{$rg}{'median'} = $newline;
	}
	if ($newline =~ /Standard deviation/) {
		my ($trash, $rg) = split("\t", $newline);
		$newline = <FILE>;
		chomp $newline;
		$is{$rg}{'sd'}  = $newline;
		$is{$rg}{'isu'} = $is{$rg}{'median'} + $is{$rg}{'sd'} * $sd_cutoff_cl;
		$is{$rg}{'isu'} = int($is{$rg}{'isu'}) + 1;
		$is{$rg}{'isd'} = $is{$rg}{'median'} - $is{$rg}{'sd'} * $sd_cutoff_cl;
		$is{$rg}{'isd'} = int($is{$rg}{'isd'}) - 1;

		if ($is{'isu'}) {
			$is{'isu'} = $is{$rg}{'isu'} if ($is{$rg}{'isu'} > $is{'isu'});
		} else {
			$is{'isu'} = $is{$rg}{'isu'};
		}

		#print "$rg\t$is{$rg}{'rl'}\t$is{$rg}{'median'}\t$is{$rg}{'sd'}\t$is{$rg}{'isu'}\t$is{$rg}{'isd'}\n";
	}
}
my $samtools_command;
if (defined($opts{S})) {
	$samtools_path .= '/' unless ($samtools_path =~ /\/$/);
	$samtools_command = $samtools_path . 'samtools';
} else {
	$samtools_command = 'samtools';
}
my $bwa_command;
if (defined($opts{W})) {
	$bwa_path .= '/' unless ($bwa_path =~ /\/$/);
	$bwa_command = $bwa_path . 'bwa';
} else {
	$bwa_command = 'bwa';
}
my $formatdb_command;
if (defined($opts{B})) {
	$blastall_path .= '/' unless ($blastall_path =~ /\/$/);
	$formatdb_command = $blastall_path . 'formatdb';
} else {
	$formatdb_command = 'formatdb';
}
my $blastall_command;
if (defined($opts{B})) {
	$blastall_path .= '/' unless ($blastall_path =~ /\/$/);
	$blastall_command = $blastall_path . 'blastall';
} else {
	$blastall_command = 'blastall';
}
$Bin =~ /scripts/;
my $bin_path          = $` . 'bin/';
my $bamreader_command = $bin_path . 'bamreader';
#my $reference_db      = Bio::DB::Fasta->new($reference_path) if ($alg or $blast_bp);

my $time0  = time;
my $local0 = localtime($time0);
print STDERR "$local0 Meerkat $version started\n";

&discord($bamfile, $prefix, $is_files_dir, $is_actual_dir, $threads_bwa, $blacklistfile, $clip, $bamreader_command, $sd_cutoff_disc, $isinfofile, $samtools_command, $remove_dup, $black_rg) if ($discordant);
my $time  = time;
my $local = localtime($time);
print STDERR "$local extracted discordant reads\n" if ($discordant);

&cluster($prefix, $is_files_dir, $is_actual_dir, $threads_bwa, $blacklistfile, \%is, $cut_sr, $ad_align, $alt_map_max, $alt_map_max_clip, $clip, $samtools_command, $cl, $bamreader_command) if ($cluster);
$time  = time;
$local = localtime($time);
print STDERR "$local called discordant clusters\n" if ($cluster);

&mpd($prefix, $is_files_dir, $is_actual_dir, $threads_bwa, $bamreader_command, $sd_cutoff_disc, $support_mps, $support_mpf) if ($mpd);
$time  = time;
$local = localtime($time);
print STDERR "$local called events by read pairs\n" if ($mpd);

&alg($bamreader_command, $prefix, $is_files_dir, $is_actual_dir, $threads_bwa, $sd_cutoff_disc, $support_mps, $support_mpf, $remove_dup, $clip, $cut_sr, $reference_path) if ($alg);
$time  = time;
$local = localtime($time);
print STDERR "$local aligned split reads\n" if ($alg);

&srd($bamreader_command, $prefix, $is_files_dir, $is_actual_dir, $threads_bwa, $sd_cutoff_disc, $support_mps, $support_mpf, $remove_dup, $clip, $reference_path, $cut_sr, $support_reads) if ($srd);
$time  = time;
$local = localtime($time);
print STDERR "$local confirmed events by split reads\n" if ($srd);

&filter($bamreader_command, $prefix, $is_files_dir, $is_actual_dir, $threads_bwa) if ($filter);
$time  = time;
$local = localtime($time);
print STDERR "$local filtered events\n" if ($filter);
&blast_bp($bamreader_command, $prefix, $is_files_dir, $is_actual_dir, $threads_bwa, $sd_cutoff_disc, $support_mps, $support_mpf, $remove_dup, $clip, $reference_path, $cut_sr, $support_reads) if ($blast_bp);
$time  = time;
$local = localtime($time);
print STDERR "$local refined break points by local alignments\n" if ($blast_bp);
my $time1      = time;
my $difference = $time1 - $time0;
my $seconds    = $difference % 60;
$difference = ($difference - $seconds) / 60;
my $minutes = $difference % 60;
$difference = ($difference - $minutes) / 60;
print STDERR "Time used: $difference:$minutes:$seconds\n";

#	blast break points reads to break points regions
	#	file format of prefix.sr.intra.refined
	#	del, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, homology sizes
	#	del_ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr (donor), range of insertion (2 col), insert size, distance of deletion and insertion
	#	del_invers, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr (donor), range of inversion (2 col), inversion size, distance of inverstion and deletion (2 col), homology at break points
	#	ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, rchr (donor), ange of insertion (2 col), insert size, distance of deletion and insertion, homology at break points
	#	invers, cluster id, number of supporting read pairs, number of supporting split reads, chr, inversion left boundary, inversion right boundary, inversion size, homology at break points
	#	invers_*, cluster id, number of supporting read pairs, number of supporting split reads, chr, inversion left boundary, inversion right boundary, inversion size, homology at break points
	#	tandem_dup, cluster id, number of supporting read pairs, number of supporting split reads, chr, tandem duplication boundary 1, tandem duplication boundary 2, tandem duplication size, homology at break points
	#	file format of prefix.sr.inter.refine
	#	del_ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr of insertion donor, range of insertion (2 col), insert size, homology at break points
	#	ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of insert site, insert site size, chr of insertion donor, range of insertion (2 col), insert size, homology at break points
	#	transl_inter, cluster id, number of supporting read pairs, number of supporting split reads, chr of 1st cluster, boundary of 1st cluster, orientation of 1st cluster, chr of 2nd cluster, boundary of 2nd cluster, orientation of 2nd cluster, homology at break points

sub blast_bp {
	my $bamreader_command = shift;
	my $prefix          = shift;
	my $working_prefix = shift;
	my $is_actual_dir = shift;
	my $threads_bwa = shift;
	my $sd_cutoff_disc   = shift;
	my $support_mps      = shift;
	my $support_mpf      = shift;
	my $remove_dup       = shift;
	my $clip             = shift;
	my $reference_path   = shift;
	my $cut_sr          = shift;
	my $support_reads   = shift;
	if("null" ne $tmp_dir) {
		print("$bamreader_command blast_bp -t $threads_bwa -x $prefix -D $is_actual_dir -d $sd_cutoff_disc -p $support_mps -o $support_mpf -m $remove_dup -l $clip -F $reference_path -s $cut_sr -q $support_reads\n");
		system("$bamreader_command blast_bp -t $threads_bwa -x $prefix -D $is_actual_dir -d $sd_cutoff_disc -p $support_mps -o $support_mpf -m $remove_dup -l $clip -F $reference_path -s $cut_sr -q $support_reads");	
	} else {
		print("$bamreader_command blast_bp -t $threads_bwa -x $prefix -d $sd_cutoff_disc -p $support_mps -o $support_mpf -m $remove_dup -l $clip -F $reference_path -s $cut_sr -q $support_reads\n");
		system("$bamreader_command blast_bp -t $threads_bwa -x $prefix -d $sd_cutoff_disc -p $support_mps -o $support_mpf -m $remove_dup -l $clip -F $reference_path -s $cut_sr -q $support_reads");
	}
}

#	filter outputs
sub filter {
	my $bamreader_command = shift;
	my $prefix         = shift;
	my $working_prefix = shift; 
	my $is_actual_dir = shift;
	my $threads_bwa = shift;
	if("null" ne $tmp_dir) {
		print("$bamreader_command filter -t $threads_bwa -x $prefix -D $is_actual_dir\n");
		system("$bamreader_command filter -t $threads_bwa -x $prefix -D $is_actual_dir");
	} else {
		print("$bamreader_command filter -t $threads_bwa -x $prefix\n");
		system("$bamreader_command filter -t $threads_bwa -x $prefix");
	}
}

#	call all SVs based on split reads
	#	file format of prefix.sr.intra.out
	#	del, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size
	#	del_ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr (donor), range of insertion (2 col), insert size, distance of deletion and insertion
	#	del_invers, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr (donor), range of inversion (2 col), inversion size, distance of inverstion and deletion (2 col)
	#	ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr (donor), range of insertion (2 col), insert size, distance of deletion and insertion
	#	invers, cluster id, number of supporting read pairs, number of supporting split reads, chr, inversion left boundary, inversion right boundary, inversion size
	#	invers_*, cluster id, number of supporting read pairs, number of supporting split reads, chr, inversion left boundary, inversion right boundary, inversion size
	#	tandem_dup, cluster id, number of supporting read pairs, number of supporting split reads, chr, tandem duplication boundary 1, tandem duplication boundary 2, tandem duplication size
	#	file format of prefix.sr.inter.out
	#	del_ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr of insertion donor, range of insertion (2 col), insert size
	#	ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of insert site, insert site size, chr of insertion donor, range of insertion (2 col), insert size
	#	transl_inter, cluster id, number of supporting read pairs, number of supporting split reads, chr of 1st cluster, boundary of 1st cluster, orientation of 1st cluster, chr of 2nd cluster, boundary of 2nd cluster, orientation of 2nd cluster
	#	$left_bound_sr1: left bound in primary cluster
	#	$right_bound_sr1: right bound in primary cluster
	#	$left_bound_sr2: left bound in secondary cluster
	#	$right_bound_sr2: right bound in secondary cluster
	#	primary cluster: cluster to initiate call in mpd
	#	secondary cluster: cluster to accompany primary cluster in mpd

sub srd {
	my $bamreader_command = shift;
	my $prefix          = shift;
	my $working_prefix = shift;
	my $is_actual_dir = shift;
	my $threads_bwa = shift;
	my $sd_cutoff_disc   = shift;
	my $support_mps      = shift;
	my $support_mpf      = shift;
	my $remove_dup       = shift;
	my $clip             = shift;
	my $reference_path   = shift;
	my $cut_sr          = shift;
	my $support_reads   = shift;
	if("null" ne $tmp_dir) {
		print("$bamreader_command srd -t $threads_bwa -x $prefix -D $is_actual_dir -d $sd_cutoff_disc -p $support_mps -o $support_mpf -m $remove_dup -l $clip -F $reference_path -s $cut_sr -q $support_reads\n");
		system("$bamreader_command srd -t $threads_bwa -x $prefix -D $is_actual_dir -d $sd_cutoff_disc -p $support_mps -o $support_mpf -m $remove_dup -l $clip -F $reference_path -s $cut_sr -q $support_reads");
	} else {
		print("$bamreader_command srd -t $threads_bwa -x $prefix -d $sd_cutoff_disc -p $support_mps -o $support_mpf -m $remove_dup -l $clip -F $reference_path -s $cut_sr -q $support_reads\n");
		system("$bamreader_command srd -t $threads_bwa -x $prefix -d $sd_cutoff_disc -p $support_mps -o $support_mpf -m $remove_dup -l $clip -F $reference_path -s $cut_sr -q $support_reads");
	}
}

#	construct search space and align split reads
sub alg {
	my $bamreader_command = shift;
	my $prefix            = shift;
	my $working_prefix = shift;
	my $is_actual_dir = shift;
	my $threads_bwa = shift;
	my $sd_cutoff_disc    = shift;
	my $support_mps       = shift;
	my $support_mpf       = shift;
	my $remove_dup        = shift;
	my $clip              = shift;
	my $cut_sr            = shift;
	my $reference_path    = shift;
	if("null" ne $tmp_dir) {
		print("$bamreader_command alg -x $prefix -D $is_actual_dir -t $threads_bwa -d $sd_cutoff_disc -p $support_mps -o $support_mpf -m $remove_dup -l $clip -s $cut_sr -F $reference_path\n");
		system("$bamreader_command alg -x $prefix -D $is_actual_dir -t $threads_bwa -d $sd_cutoff_disc -p $support_mps -o $support_mpf -m $remove_dup -l $clip -s $cut_sr -F $reference_path");
	} else {
		print("$bamreader_command alg -x $prefix -t $threads_bwa -d $sd_cutoff_disc -p $support_mps -o $support_mpf -m $remove_dup -l $clip -s $cut_sr -F $reference_path\n");
		system("$bamreader_command alg -x $prefix -t $threads_bwa -d $sd_cutoff_disc -p $support_mps -o $support_mpf -m $remove_dup -l $clip -s $cut_sr -F $reference_path");
	}
}

#	extract discordant read pairs
sub discord {
	my $bamfile          = shift;
	my $prefix           = shift;
	my $working_prefix = shift;
	my $is_actual_dir = shift;
	my $threads_bwa = shift;
	my $blacklistfile    = shift;
	my $clip             = shift;
	my $bamreader_command = shift;
	my $sd_cutoff_disc   = shift;
	my $isinfofile       = shift;
	my $samtools_command = shift;
	my $remove_dup       = shift;
	my $black_rg         = shift;
	my $drelog           = $prefix . '.dre.log';
	my $cl_bam           = $prefix . '.cl.sorted.bam';
	my $disc_bam         = $prefix . '.disc.bam';
	my $dup_file         = $prefix . '.dup.bam';
	my $disc_cl_bam      = $prefix . '.cl.disc.bam';
	my $dup_cl_file      = $prefix . '.cl.dup.bam';
	my $disc_sort    = $prefix . '.disc.sorted.bam';
	my $disc_sort_num    = $prefix . '.disc.sorted.num.bam';
	my $disc_cl_sort = $prefix . '.cl.disc.sorted.bam';
	my $disc_cl_sort_num = $prefix . '.cl.disc.sorted.num.bam';
	if("null" ne $tmp_dir) {
		if(-e $working_prefix) {
			$drelog           = $working_prefix. '.dre.log';
			$cl_bam           = $working_prefix. '.cl.sorted.bam';
			$disc_bam         = $working_prefix. '.disc.bam';
			$dup_file         = $working_prefix. '.dup.bam';
			$disc_cl_bam      = $working_prefix. '.cl.disc.bam';
			$dup_cl_file      = $working_prefix. '.cl.dup.bam';
			$disc_sort    = $working_prefix. '.disc.sorted.bam';
			$disc_sort_num    = $working_prefix. '.disc.sorted.num.bam';
			$disc_cl_sort = $working_prefix. '.cl.disc.sorted.bam';
			$disc_cl_sort_num = $working_prefix. '.cl.disc.sorted.num.bam';
		}
	}
	die "bam file doesn't exist\n" unless (-e $bamfile);

	if ($blacklistfile) {
		if ($remove_dup) {
			if("null" ne $tmp_dir) {
				print("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -P -s $sd_cutoff_disc  -r $black_rg -i $isinfofile -o $disc_bam -d $dup_file -a $bamfile >$drelog\n");
				system("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -P -s $sd_cutoff_disc  -r $black_rg -i $isinfofile -o $disc_bam -d $dup_file -a $bamfile >$drelog");
			} else {
				print("$bamreader_command dre -t $threads_bwa -f -v -P -s $sd_cutoff_disc  -r $black_rg -i $isinfofile -o $disc_bam -d $dup_file -a $bamfile >$drelog\n");
				system("$bamreader_command dre -t $threads_bwa -f -v -P -s $sd_cutoff_disc  -r $black_rg -i $isinfofile -o $disc_bam -d $dup_file -a $bamfile >$drelog");
			}
			my $time  = time;
			my $local = localtime($time);
			print STDERR "$local extracted discordant read pairs\n";
			if ($clip) {
				# some times this blacklist file causes the segment fault error since it is disabled for small .cl.bam files.
				if("null" ne $tmp_dir) {
#					print("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -P -s $sd_cutoff_disc  -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog\n");
#					system("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -P -s $sd_cutoff_disc  -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog");
					print("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog\n");
					system("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog");
				} else {
#					print("$bamreader_command dre -t $threads_bwa -f -v -P -s $sd_cutoff_disc  -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog\n");
#					system("$bamreader_command dre -t $threads_bwa -f -v -P -s $sd_cutoff_disc  -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog");
					print("$bamreader_command dre -t $threads_bwa -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog\n");
					system("$bamreader_command dre -t $threads_bwa -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog");
				}
				my $time  = time;
				my $local = localtime($time);
				print STDERR "$local extracted discordant read pairs from remapped bam file\n";
			}
		} else {
			if("null" ne $tmp_dir) {
				# 
				print("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -m -s $sd_cutoff_disc  -r $black_rg -i $isinfofile -o $disc_bam -a $bamfile >$drelog\n");
				system("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -m -s $sd_cutoff_disc  -r $black_rg -i $isinfofile -o $disc_bam -a $bamfile >$drelog");
			} else {
				print("$bamreader_command dre -t $threads_bwa -f -v -m -s $sd_cutoff_disc  -r $black_rg -i $isinfofile -o $disc_bam -a $bamfile >$drelog\n");
				system("$bamreader_command dre -t $threads_bwa -f -v -m -s $sd_cutoff_disc  -r $black_rg -i $isinfofile -o $disc_bam -a $bamfile >$drelog");
			}
			my $time  = time;
			my $local = localtime($time);
			print STDERR "$local extracted discordant read pairs\n";
			if ($clip) {
				# some times this blacklist file causes the segment fault error since it is disabled for small .cl.bam files.
				if("null" ne $tmp_dir) {
#					print("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -P -s $sd_cutoff_disc  -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog\n");
#					system("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -P -s $sd_cutoff_disc  -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog");
					print("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog\n");
					system("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog");
				} else {
#					print("$bamreader_command dre -t $threads_bwa -f -v -P -s $sd_cutoff_disc  -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog\n");
#					system("$bamreader_command dre -t $threads_bwa -f -v -P -s $sd_cutoff_disc  -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog");
					print("$bamreader_command dre -t $threads_bwa -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog\n");
					system("$bamreader_command dre -t $threads_bwa -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog");
				}
				my $time  = time;
				my $local = localtime($time);
				print STDERR "$local extracted discordant read pairs from remapped bam file\n";
			}
		}
	} else {
		if ($remove_dup) {
			if("null" ne $tmp_dir) {
				print("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_bam -d $dup_file -a $bamfile >$drelog\n");
				system("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_bam -d $dup_file -a $bamfile >$drelog");
			} else {
				print("$bamreader_command dre -t $threads_bwa -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_bam -d $dup_file -a $bamfile >$drelog\n");
				system("$bamreader_command dre -t $threads_bwa -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_bam -d $dup_file -a $bamfile >$drelog");
			}
			my $time  = time;
			my $local = localtime($time);
			print STDERR "$local extracted discordant read pairs\n";
			if ($clip) {
				if("null" ne $tmp_dir) {
					print("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog\n");
					system("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog");
				} else {
					print("$bamreader_command dre -t $threads_bwa -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog\n");
					system("$bamreader_command dre -t $threads_bwa -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog");
				}
				my $time  = time;
				my $local = localtime($time);
				print STDERR "$local extracted discordant read pairs from remapped bam file\n";
			}
		} else {
			if("null" ne $tmp_dir) {
				print("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -m -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_bam -a $bamfile >$drelog\n");
				system("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -m -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_bam -a $bamfile >$drelog");
			} else {
				print("$bamreader_command dre -t $threads_bwa -f -v -m -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_bam -a $bamfile >$drelog\n");
				system("$bamreader_command dre -t $threads_bwa -f -v -m -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_bam -a $bamfile >$drelog");
			}
			my $time  = time;
			my $local = localtime($time);
			print STDERR "$local extracted discordant read pairs\n";
			if ($clip) {
				if("null" ne $tmp_dir) {
					print("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog\n");
					system("$bamreader_command dre -t $threads_bwa -D $is_actual_dir -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog");
				} else {
					print("$bamreader_command dre -t $threads_bwa -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog\n");
					system("$bamreader_command dre -t $threads_bwa -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file -a $cl_bam >>$drelog");
				}
				my $time  = time;
				my $local = localtime($time);
				print STDERR "$local extracted discordant read pairs from remapped bam file\n";
			}
		}
	}
#	system "$samtools_command sort -n -@ $threads_bwa $disc_bam $disc_sort";
#	print("sambamba sort -N -l 1 -t $threads_bwa -o $disc_sort $disc_bam\n");
#	system("sambamba sort -N -l 1 -t $threads_bwa -o $disc_sort $disc_bam");
#	print("sambamba sort -l 1 -t $threads_bwa -o $disc_sort_num $disc_bam\n");
#	system("sambamba sort -l 1 -t $threads_bwa -o $disc_sort_num $disc_bam");
	
	#	system "rm $disc_bam";
	#	my $time  = time;
	#	my $local = localtime($time);
	#	print STDERR "$local sorted discordant bam by name\n";
#	if ($clip) {
#		system "$samtools_command sort -n -@ $threads_bwa $disc_cl_bam $disc_cl_sort";
#		print("sambamba sort -N -l 1 -t $threads_bwa -o $disc_cl_sort $disc_cl_bam\n");
#		system("sambamba sort -N -l 1 -t $threads_bwa -o $disc_cl_sort $disc_cl_bam");
#		print("sambamba sort -l 1 -t $threads_bwa -o $disc_cl_sort_num $disc_cl_bam\n");
#		system("sambamba sort -l 1 -t $threads_bwa -o $disc_cl_sort_num $disc_cl_bam");
		#		system "rm $disc_cl_bam";
#		my $time  = time;
#		my $local = localtime($time);
#		print STDERR "$local sorted discordant bam by name from remapped bam file\n";
#	}
}

#	cluster discordant read pairs
	#	prefix.discord file format: cluster id, number of read pairs, chr, position, strand, mchr, mposition, mstrand, size
	#	prefix.clusters file format: cluster id primary, cluster id secondary, old cluster id, weight, mismatch, readname, read group, chr, strand, start, length, mchr, mstrand, mstart, mlength, isize
	# hash keys
	# n: readname
	# c: seqid
	# t: start
	# e: end
	# d: strand
	# m: nm
	# l: len
	# C: mseqid
	# T: mstart
	# E: mend
	# D: mstrad
	# M: mnm
	# L: mlen
	# i: isize
	# r: rg
	# w: weight
	# a: alternative mapping

sub cluster {
	my $prefix              = shift;
	my $working_prefix = shift;
	my $is_actual_dir = shift;
	my $threads_bwa = shift;
	my $blacklistfile       = shift;
	my $ref_is              = shift;
	my $cut_sr              = shift;
	my $ad_align            = shift;
	my $alt_map_max         = shift;
	my $alt_map_max_clip    = shift;
	my $clip                = shift;
	my $samtools_command    = shift;
	my $cl                  = shift;
	my $bamreader_command    = shift;
	my $disc_bam            = $prefix . '.disc.sorted.bam';
	my $disc_cl_bam         = $prefix . '.cl.disc.sorted.bam';
	my $raw_alt_map_file    = $prefix . '.alt.map.raw';
	my $raw_alt_map_cl_file = $prefix . '.alt.map.raw.cl';
	my $sort_alt_map_file   = $prefix . '.alt.map.sort';
	my $raw_mapping         = $prefix . '.mapping.raw';
	my $raw_cluster         = $prefix . '.clusters.raw';         # for -a 0 only
	my $clusterfile         = $prefix . '.clusters';             # final output
	if("null" ne $tmp_dir) {
		if(-e $working_prefix) {
			$disc_bam            = $working_prefix. '.disc.sorted.bam';
			$disc_cl_bam         = $working_prefix. '.cl.disc.sorted.bam';
			$raw_alt_map_file    = $working_prefix. '.alt.map.raw';
			$raw_alt_map_cl_file = $working_prefix. '.alt.map.raw.cl';
			$sort_alt_map_file   = $working_prefix. '.alt.map.sort';
			$raw_mapping         = $working_prefix. '.mapping.raw';
			$raw_cluster         = $working_prefix. '.clusters.raw';         # for -a 0 only
			$clusterfile         = $working_prefix. '.clusters';             # final output
		}
	}
#	goto AS2 if ($cl == 2 or $cl == 3);

	# parse positions to be ignored
	my %blacklist;
	if (-e $blacklistfile) {
		open(FILE, "gunzip -c $blacklistfile |");
		my $newline;
		while ($newline = <FILE>) {
			chomp $newline;
			my @data = split(/\t/, $newline);
			$blacklist{$data[0]}{$data[1]} = 1;

			#print "$data[0]\t$data[1]\n";
		}
		close FILE;
	}
	if("null" ne $tmp_dir) {
		print("$bamreader_command r_alt -t $threads_bwa -x $prefix -D $is_actual_dir -f $alt_map_max -s $cut_sr -i $disc_bam -z $raw_alt_map_file\n");
		system("$bamreader_command r_alt -t $threads_bwa -x $prefix -D $is_actual_dir -f $alt_map_max -s $cut_sr -i $disc_bam -z $raw_alt_map_file");
	} else {
		print("$bamreader_command r_alt -t $threads_bwa -x $prefix -f $alt_map_max -s $cut_sr -i $disc_bam -z $raw_alt_map_file\n");
		system("$bamreader_command r_alt -t $threads_bwa -x $prefix -f $alt_map_max -s $cut_sr -i $disc_bam -z $raw_alt_map_file");
	}
#	if ($clip) {
#		if("null" ne $tmp_dir) {
#			system("$bamreader_command r_alt -t $threads_bwa -x $prefix -D $is_actual_dir -f $alt_map_max_clip -s $cut_sr -i $disc_cl_bam -z $raw_alt_map_cl_file\n");
#			system("$bamreader_command r_alt -t $threads_bwa -x $prefix -D $is_actual_dir -f $alt_map_max_clip -s $cut_sr -i $disc_cl_bam -z $raw_alt_map_cl_file");
#		} else {
#			system("$bamreader_command r_alt -t $threads_bwa -x $prefix -f $alt_map_max_clip -s $cut_sr -i $disc_cl_bam -z $raw_alt_map_cl_file\n");
#			system("$bamreader_command r_alt -t $threads_bwa -x $prefix -f $alt_map_max_clip -s $cut_sr -i $disc_cl_bam -z $raw_alt_map_cl_file");
#		}
#		system "cat $raw_alt_map_cl_file >> $raw_alt_map_file";
##		system "rm $raw_alt_map_cl_file";
#	}

	#	close $rawalt_fh;
	my $time  = time;
	my $local = localtime($time);
	print STDERR "$local constructed all mappings\n";
	print("sort -k 2,2 -k 3,3n $raw_alt_map_file > $sort_alt_map_file\n");
	system("sort -k 2,2 -k 3,3n $raw_alt_map_file > $sort_alt_map_file");

	#system "rm $raw_alt_map_file";
#	return if ($cl == 1);

	#	construct clusters for all alternative mappings include uniq mapped reads
#  AS2: goto AS3 if ($cl == 3);

	#	%cluster_alt: cluster of all alternative mappings, same format as %cluster
	#system "rm $sort_alt_map_file";
	if("null" ne $tmp_dir) {
		print("$bamreader_command s_alt -t $threads_bwa -x $prefix -D $is_actual_dir -c $sd_cutoff_cl $bamfile\n");
		system("$bamreader_command s_alt -t $threads_bwa -x $prefix -D $is_actual_dir -c $sd_cutoff_cl $bamfile");
	} else {
		print("$bamreader_command s_alt -t $threads_bwa -x $prefix -c $sd_cutoff_cl $bamfile\n");
		system("$bamreader_command s_alt -t $threads_bwa -x $prefix -c $sd_cutoff_cl $bamfile");
	}
	$time  = time;
	$local = localtime($time);
	print STDERR "$local constructed discordant clusters for all alternative mappings\n";
#	return if ($cl == 2);
#  AS3: 
  if ($ad_align) {
  	print("$bamreader_command sclus -i $raw_mapping -z $clusterfile\n");
	system("$bamreader_command sclus -t $threads_bwa -i $raw_mapping -z $clusterfile");

		#system "rm $raw_mapping";
		#		my $time  = time;
		#		my $local = localtime($time);
		#		print STDERR "$local called discordant clusters\n";
	} else {
		system "sort -k 1,1n $raw_mapping > $raw_cluster";
		open CLFILE, ">$clusterfile";
		open RAWCL,  "<$raw_cluster";
		while ($newline = <RAWCL>) {
			chomp $newline;
			my @data = split(/\t/, $newline);
			print CLFILE "$data[0]\t0\t\t\t\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$data[5]\t$data[6]\t$data[7]\t$data[8]\t$data[9]\t$data[10]\t$data[11]\n";
		}
		close CLFILE;

		#system "rm $raw_cluster";
		#system "rm $raw_mapping";
		#		my $time  = time;
		#		my $local = localtime($time);
		#		print STDERR "$local called discordant clusters\n";
	}
}

#	identify candidates from read pairs

	#	file format of prefix.mp.intra.out
	#	del: simple deletion
	#	del_inssd: deletion with insertion in the same orientation, insertion come from downstream of deletion
	#	del_inssu: deletion with insertion in the same orientation, insertion come from upstream of deletion
	#	del_insod: deletion with insertion in the opposite orientation, insertion come from downstream of deletion
	#	del_insou: deletion with insertion in the opposite orientation, insertion come from upstream of deletion
	#	invers: inversion supported by 2 reciprocal clusters, include simple inversion and deletion with inversion
	#	tandem_dup: tandem duplication
	#	invers_f: inversion supported by one cluster on forward strand
	#	invers_r: inversion supported by one cluster on reverse strand
	#	del, cluster id, number of supporting read pairs, chr, range of deletion (2 col), deletion size, region of break point
	#	del_ins*, cluster id, number of supporting read pairs, chr, range of deletion (2 col), deletion size, chr (donor), range of insertion (2 col), insert size, distance of deletion and insertion, region of break point (2 col)
	#	*invers, cluster id, number of supporting read pairs, chr, inversion left boundary (2 col), inversion size, inversion right boundary (2 col), region of break point (2 col)
	#	invers_*, cluster id, number of supporting read pairs, chr, inversion left boundary, inversion right boundary, inversion size, region of break point
	#	tandem_dup, cluster id, number of supporting read pairs, chr, tandem duplication boundary 1, tandem duplication boundary 2, tandem duplication size, region of break point
	#	file format of prefix.mp.inter.out
	#	del_inss: deletion with insertion in the same orientation
	#	del_inso: deletion with insertion in the opposite orientation
	#	transl_inter: inter-chromosomal translocation
	#	del_ins*, cluster id, number of supporting read pairs, chr of deletion, range of deletion (2 col), deletion size, chr of insertion donor, range of insertion (2 col), insert size, region of break point (2 col)
	#	transl_inter, cluster id, number of supporting read pairs, chr of 1st cluster, range (2 col), orientation of 1st cluster, chr of 2nd cluster, range (2 col), orientation of 2nd cluster, region of break point


sub mpd {
	my $prefix         = shift;
	my $working_prefix = shift;
	my $is_actual_dir = shift;
	my $threads_bwa = shift;
	my $bamreader      = shift;
	my $sd_cutoff_disc = shift;
	my $support_mps    = shift;
	my $support_mpf    = shift;
	if("null" ne $tmp_dir) {
		print("$bamreader mpd -t $threads_bwa -x $prefix -D $is_actual_dir -d $sd_cutoff_disc -p $support_mps -o $support_mpf\n");
		system("$bamreader mpd -t $threads_bwa -x $prefix -D $is_actual_dir -d $sd_cutoff_disc -p $support_mps -o $support_mpf");
	} else {
		print("$bamreader mpd -t $threads_bwa -x $prefix -d $sd_cutoff_disc -p $support_mps -o $support_mpf\n");
		system("$bamreader mpd -t $threads_bwa -x $prefix -d $sd_cutoff_disc -p $support_mps -o $support_mpf");
	}
}

sub print_usage {
	die "Usage:\nperl ./scripts/meerkat.pl [options]
	-b FILE	sorted and indexed bam file, required
	-k INT	[0/1], use black list generated in pre_process.pl, default 1
	-d FLT	standard deviation cutoff to call discordant read pairs, default 3
	-c FLT	standard deviation cutoff to cluster discordant read pairs, default equal to -d, it is recommended to use same -c as -d if -d<=5, -c 5 if -d > 5
	-p INT	number of supporting read pairs required for an event to be called, default 2
	-o INT	number of supporting full length read pairs, default 0, specify this option will decrease sensitivity on small complex events
	-q INT	number of supporting split reads required, default 1
	-z INT	event size cutoff, default 1,000,000,000
	-s INT	bp to be cut off from beginning and end of unmapped reads, must be same as -s in pre_process.pl, default 20
	-m INT	[0/1], if set to 1, use Meerkat to remove duplicates; if set to 0, use flag 'd' marked by Picard or other tools to remove duplicates. If the bam file is aligned by bwa mem, it has to be processed by Picard to mark duplicates and use 0 option. The bwa mem aligned bam file won't work with option 1. Default 1
	-a INT	[0/1], adjust non-uniq mapped reads, default 1
	-u INT  [0/1], use all alignments in the BAM file, turn this option on if the BAM file is not generated by BWA, turn on this option will force turning off option a, default 0
	-Q INT	minimum mapping quality for reads to be used, default 0
	-g INT	number of alternative mappings to consider in main bam file, number of alternative mappings printed out in XA tag by bwa is controlled by -N, default use all in bam file
	-f INT	number of alternative mappings to consider in clipped alignments, default use all in bam file
	-l INT	[0/1], consider clipped alignments, default 1
	-t INT	number of threads used in bwa alignment, default 1
	-R STR	file name of read group to be ignored, one read group ID per line
	-F STR	/path/to/reference/fasta/files, path only, not the files, required
	-S STR	/path/to/samtools, path only, not the command, no need to specify if samtools is in PATH
	-W STR	/path/to/bwa, path only, not the command, no need to specify if bwa is in PATH
	-B STR	/path/to/blastall and formatdb, path only, not the command, no need to specify if blastall and formatdb is in PATH
	-P STR	specify step to run, dc|cl|mpd|alg|srd|rf|all, default 'all', each step require results from previous steps
		dc: extract discordant read pairs
		cl: construct clusters of discordant read pairs
		mpd: call events based on read pairs
		alg: align split reads to candidate break point regions
		srd: confirm events based on split reads and filter results
		rf:  refine break points by local alignments
		all: run all above steps
	-h help\n";
}

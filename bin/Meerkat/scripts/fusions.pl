#!/usr/bin/perl
# fusions.pl
# identify gene that has been mutated with given variants produced by Meerkat and predict fusion types
# Author: Lixing Yang, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, USA
# Email: lixing_yang@hms.harvard.edu
# i.e. perl scripts/fusions.pl -i inputfile -G gene_annotation


use strict;
use Getopt::Std;

my %opts = ();
getopts("i:G:h", \%opts);

my $inputfile = $opts{i};
my $refseq_gene = $opts{G};
if (defined($opts{h}))
{
	&print_usage;
	die;
}

die "Please provide gene annotation file\nFor example download human reference gene annotation from \nhttp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz\nand unzip it.\n" unless (defined($refseq_gene));

my $newline;
my (%gene, %exon, %exon2, %count);
my $i = 0;
open FILE, "<$refseq_gene";
while ($newline = <FILE>)
{
	chomp $newline;
	my @data = split (/\t/, $newline);
	$gene{$data[2]}[$i][0] = $data[12]; #gene name
	$gene{$data[2]}[$i][1] = $data[4]; #gene start
	$gene{$data[2]}[$i][2] = $data[5]; #gene end
	$i++;
	$exon{$data[12]}[0] = $data[3]; #orientation
	$count{$data[12]} = 0 unless (defined($count{$data[12]}));
	my @start = split (/,/, $data[9]);
	my @end = split (/,/, $data[10]);
	my @frame = split (/,/, $data[15]);
	my $j = 0;
	foreach (@start)
	{
		my @tmp = ($start[$j], $end[$j], $frame[$j]);
		push @{$exon{$data[12]}[1]}, \@tmp;
		push @{$exon2{$data[12]}[$count{$data[12]}][1]}, \@tmp; # codon frame
		$j++;
	}
	$exon2{$data[12]}[$count{$data[12]}][0] = $data[3]; #orientation
	$exon2{$data[12]}[$count{$data[12]}][2] = $data[4]; #gene start
	$exon2{$data[12]}[$count{$data[12]}][3] = $data[5]; #gene end
	$exon2{$data[12]}[$count{$data[12]}][4] = $data[6]-1; #cds start
	$exon2{$data[12]}[$count{$data[12]}][5] = $data[7]+1; #cds end
	$count{$data[12]}++;
}
close FILE;

die "Input must be *.variants file generated by Meerkat\n" unless ($inputfile =~ /.variants/);
my $mutgenefile = $`.'.mut';
my $fusiontypefile = $`.'.fusions';
my @variant;
open FILE, "<$inputfile";
while ($newline = <FILE>)
{
	chomp $newline;
	my @data = split (/\t/, $newline);
	push @variant, \@data;
}
close FILE;
	
open OUT, ">$mutgenefile";
foreach (@variant)
{
	my $variant = $_;
	my $toprint;
	$$variant[-1] =~ /RP:/;
	my @rp = split (/\//, $');#'
	#print "@rp\n";
	# intra
	if ($$variant[0] eq 'del')
	{
		my ($gene1, $exon_intron1, $up_gene1, $up_distance1, $down_gene1, $down_distance1) = &mut_gene($$variant[5], $$variant[6]);
		my ($gene2, $exon_intron2, $up_gene2, $up_distance2, $down_gene2, $down_distance2) = &mut_gene($$variant[5], $$variant[7]);
		$toprint = "$$variant[5]\t$$variant[6]\t1\t$gene1\t$exon_intron1\t$$variant[5]\t$$variant[7]\t-1\t$gene2\t$exon_intron2\t$$variant[0]\t$$variant[1]\t$$variant[2]\t$$variant[3]\t$$variant[4]\t$$variant[9]\t$rp[0]\n";
	}
	elsif ($$variant[0] eq 'del_ins')
	{
		my ($gene1, $exon_intron1, $up_gene1, $up_distance1, $down_gene1, $down_distance1) = &mut_gene($$variant[5], $$variant[6]);
		my ($gene2, $exon_intron2, $up_gene2, $up_distance2, $down_gene2, $down_distance2) = &mut_gene($$variant[5], $$variant[7]);
		$toprint = "$$variant[5]\t$$variant[6]\t1\t$gene1\t$exon_intron1\t$$variant[5]\t$$variant[7]\t-1\t$gene2\t$exon_intron2\t$$variant[0]\t$$variant[1]\t$$variant[2]\t$$variant[3]\t$$variant[4]\t-$$variant[12]\t$rp[0]\n";
	}
	elsif ($$variant[0] eq 'tandem_dup')
	{
		my ($gene1, $exon_intron1, $up_gene1, $up_distance1, $down_gene1, $down_distance1) = &mut_gene($$variant[5], $$variant[6]);
		my ($gene2, $exon_intron2, $up_gene2, $up_distance2, $down_gene2, $down_distance2) = &mut_gene($$variant[5], $$variant[7]);
		$toprint = "$$variant[5]\t$$variant[6]\t-1\t$gene1\t$exon_intron1\t$$variant[5]\t$$variant[7]\t1\t$gene2\t$exon_intron2\t$$variant[0]\t$$variant[1]\t$$variant[2]\t$$variant[3]\t$$variant[4]\t$$variant[9]\t$rp[0]\n";
	}
	elsif ($$variant[0] eq 'invers_f')
	{
		my ($gene1, $exon_intron1, $up_gene1, $up_distance1, $down_gene1, $down_distance1) = &mut_gene($$variant[5], $$variant[6]);
		my ($gene2, $exon_intron2, $up_gene2, $up_distance2, $down_gene2, $down_distance2) = &mut_gene($$variant[5], $$variant[7]);
		$toprint = "$$variant[5]\t$$variant[6]\t1\t$gene1\t$exon_intron1\t$$variant[5]\t$$variant[7]\t1\t$gene2\t$exon_intron2\t$$variant[0]\t$$variant[1]\t$$variant[2]\t$$variant[3]\t$$variant[4]\t$$variant[9]\t$rp[0]\n";
	}
	elsif ($$variant[0] eq 'invers_r')
	{
		my ($gene1, $exon_intron1, $up_gene1, $up_distance1, $down_gene1, $down_distance1) = &mut_gene($$variant[5], $$variant[6]);
		my ($gene2, $exon_intron2, $up_gene2, $up_distance2, $down_gene2, $down_distance2) = &mut_gene($$variant[5], $$variant[7]);
		$toprint = "$$variant[5]\t$$variant[6]\t-1\t$gene1\t$exon_intron1\t$$variant[5]\t$$variant[7]\t-1\t$gene2\t$exon_intron2\t$$variant[0]\t$$variant[1]\t$$variant[2]\t$$variant[3]\t$$variant[4]\t$$variant[9]\t$rp[0]\n";
	}
	elsif ($$variant[0] eq 'invers')
	{
		my ($gene1, $exon_intron1, $up_gene1, $up_distance1, $down_gene1, $down_distance1) = &mut_gene($$variant[5], $$variant[6]);
		my ($gene2, $exon_intron2, $up_gene2, $up_distance2, $down_gene2, $down_distance2) = &mut_gene($$variant[5], $$variant[7]);
		my @id = split ("\/", $$variant[2]); $id[1] = $id[0] unless (defined($id[1]));
		my @mpd = split ("\/", $$variant[3]);$mpd[1] = $mpd[0] unless (defined($mpd[1]));
		my @srd = split ("\/", $$variant[4]);$srd[1] = $srd[0] unless (defined($srd[1]));
		my @hm = split ("\/", $$variant[9]);$hm[1] = $hm[0] unless (defined($hm[1]));
		$toprint = "$$variant[5]\t$$variant[6]\t1\t$gene1\t$exon_intron1\t$$variant[5]\t$$variant[7]\t1\t$gene2\t$exon_intron2\t$$variant[0]\t$$variant[1]\t$id[0]\t$mpd[0]\t$srd[0]\t$hm[0]\t$rp[0]\n";
		$toprint .= "$$variant[5]\t$$variant[6]\t-1\t$gene1\t$exon_intron1\t$$variant[5]\t$$variant[7]\t-1\t$gene2\t$exon_intron2\t$$variant[0]\t$$variant[1]\t$id[1]\t$mpd[1]\t$srd[1]\t$hm[1]\t$rp[0]\n";
	}
	elsif ($$variant[0] =~ /inssd/)
	{
		my ($gene1, $exon_intron1, $up_gene1, $up_distance1, $down_gene1, $down_distance1) = &mut_gene($$variant[5], $$variant[7]);
		my ($gene2, $exon_intron2, $up_gene2, $up_distance2, $down_gene2, $down_distance2) = &mut_gene($$variant[9], $$variant[11]);
		my ($gene3, $exon_intron3, $up_gene3, $up_distance3, $down_gene3, $down_distance3) = &mut_gene($$variant[5], $$variant[6]);
		my ($gene4, $exon_intron4, $up_gene4, $up_distance4, $down_gene4, $down_distance4) = &mut_gene($$variant[9], $$variant[10]);
		my @id = split ("\/", $$variant[2]);
		my @mpd = split ("\/", $$variant[3]);
		my @srd = split ("\/", $$variant[4]);
		my @hm = split ("\/", $$variant[14]);
		$toprint = "$$variant[5]\t$$variant[7]\t-1\t$gene1\t$exon_intron1\t$$variant[9]\t$$variant[11]\t1\t$gene2\t$exon_intron2\t$$variant[0]\t$$variant[1]\t$id[0]\t$mpd[0]\t$srd[0]\t$hm[0]\t$rp[0]\n";
		$toprint .= "$$variant[5]\t$$variant[6]\t1\t$gene3\t$exon_intron3\t$$variant[9]\t$$variant[10]\t-1\t$gene4\t$exon_intron4\t$$variant[0]\t$$variant[1]\t$id[1]\t$mpd[1]\t$srd[1]\t$hm[1]\t$rp[1]\n";
	}
	elsif ($$variant[0] =~ /inssu/)
	{
		my ($gene1, $exon_intron1, $up_gene1, $up_distance1, $down_gene1, $down_distance1) = &mut_gene($$variant[9], $$variant[10]);
		my ($gene2, $exon_intron2, $up_gene2, $up_distance2, $down_gene2, $down_distance2) = &mut_gene($$variant[5], $$variant[6]);
		my ($gene3, $exon_intron3, $up_gene3, $up_distance3, $down_gene3, $down_distance3) = &mut_gene($$variant[9], $$variant[11]);
		my ($gene4, $exon_intron4, $up_gene4, $up_distance4, $down_gene4, $down_distance4) = &mut_gene($$variant[5], $$variant[7]);
		my @id = split ("\/", $$variant[2]);
		my @mpd = split ("\/", $$variant[3]);
		my @srd = split ("\/", $$variant[4]);
		my @hm = split ("\/", $$variant[14]);
		$toprint = "$$variant[9]\t$$variant[10]\t-1\t$gene1\t$exon_intron1\t$$variant[5]\t$$variant[6]\t1\t$gene2\t$exon_intron2\t$$variant[0]\t$$variant[1]\t$id[0]\t$mpd[0]\t$srd[0]\t$hm[0]\t$rp[0]\n";
		$toprint .= "$$variant[9]\t$$variant[11]\t1\t$gene3\t$exon_intron3\t$$variant[5]\t$$variant[7]\t-1\t$gene4\t$exon_intron4\t$$variant[0]\t$$variant[1]\t$id[1]\t$mpd[1]\t$srd[1]\t$hm[1]\t$rp[1]\n";
	}
	elsif ($$variant[0] =~ /insod/)
	{
		my ($gene1, $exon_intron1, $up_gene1, $up_distance1, $down_gene1, $down_distance1) = &mut_gene($$variant[5], $$variant[6]);
		my ($gene2, $exon_intron2, $up_gene2, $up_distance2, $down_gene2, $down_distance2) = &mut_gene($$variant[9], $$variant[11]);
		my ($gene3, $exon_intron3, $up_gene3, $up_distance3, $down_gene3, $down_distance3) = &mut_gene($$variant[5], $$variant[7]);
		my ($gene4, $exon_intron4, $up_gene4, $up_distance4, $down_gene4, $down_distance4) = &mut_gene($$variant[9], $$variant[10]);
		my @id = split ("\/", $$variant[2]);
		my @mpd = split ("\/", $$variant[3]);
		my @srd = split ("\/", $$variant[4]);
		my @hm = split ("\/", $$variant[14]);
		$toprint = "$$variant[5]\t$$variant[6]\t1\t$gene1\t$exon_intron1\t$$variant[9]\t$$variant[11]\t1\t$gene2\t$exon_intron2\t$$variant[0]\t$$variant[1]\t$id[0]\t$mpd[0]\t$srd[0]\t$hm[0]\t$rp[0]\n";
		$toprint .= "$$variant[5]\t$$variant[7]\t-1\t$gene3\t$exon_intron3\t$$variant[9]\t$$variant[10]\t-1\t$gene4\t$exon_intron4\t$$variant[0]\t$$variant[1]\t$id[1]\t$mpd[1]\t$srd[1]\t$hm[1]\t$rp[1]\n";
	}
	elsif ($$variant[0] =~ /insou/)
	{
		my ($gene1, $exon_intron1, $up_gene1, $up_distance1, $down_gene1, $down_distance1) = &mut_gene($$variant[9], $$variant[11]);
		my ($gene2, $exon_intron2, $up_gene2, $up_distance2, $down_gene2, $down_distance2) = &mut_gene($$variant[5], $$variant[6]);
		my ($gene3, $exon_intron3, $up_gene3, $up_distance3, $down_gene3, $down_distance3) = &mut_gene($$variant[9], $$variant[10]);
		my ($gene4, $exon_intron4, $up_gene4, $up_distance4, $down_gene4, $down_distance4) = &mut_gene($$variant[5], $$variant[7]);
		my @id = split ("\/", $$variant[2]);
		my @mpd = split ("\/", $$variant[3]);
		my @srd = split ("\/", $$variant[4]);
		my @hm = split ("\/", $$variant[14]);
		$toprint = "$$variant[9]\t$$variant[11]\t1\t$gene1\t$exon_intron1\t$$variant[5]\t$$variant[6]\t1\t$gene2\t$exon_intron2\t$$variant[0]\t$$variant[1]\t$id[0]\t$mpd[0]\t$srd[0]\t$hm[0]\t$rp[0]\n";
		$toprint .= "$$variant[9]\t$$variant[10]\t-1\t$gene3\t$exon_intron3\t$$variant[5]\t$$variant[7]\t-1\t$gene4\t$exon_intron4\t$$variant[0]\t$$variant[1]\t$id[1]\t$mpd[1]\t$srd[1]\t$hm[1]\t$rp[1]\n";
	}
	elsif ($$variant[0] eq 'del_invers')
	{
		my ($gene1, $exon_intron1, $up_gene1, $up_distance1, $down_gene1, $down_distance1) = &mut_gene($$variant[5], $$variant[6]);
		my ($gene2, $exon_intron2, $up_gene2, $up_distance2, $down_gene2, $down_distance2) = &mut_gene($$variant[9], $$variant[11]);
		my ($gene3, $exon_intron3, $up_gene3, $up_distance3, $down_gene3, $down_distance3) = &mut_gene($$variant[5], $$variant[7]);
		my ($gene4, $exon_intron4, $up_gene4, $up_distance4, $down_gene4, $down_distance4) = &mut_gene($$variant[9], $$variant[10]);
		my @id = split ("\/", $$variant[2]);
		my @mpd = split ("\/", $$variant[3]);
		my @srd = split ("\/", $$variant[4]);
		my @hm = split ("\/", $$variant[15]);
		$toprint = "$$variant[5]\t$$variant[6]\t1\t$gene1\t$exon_intron1\t$$variant[9]\t$$variant[11]\t1\t$gene2\t$exon_intron2\t$$variant[0]\t$$variant[1]\t$id[0]\t$mpd[0]\t$srd[0]\t$hm[0]\t$rp[0]\n";
		$toprint .= "$$variant[5]\t$$variant[7]\t-1\t$gene3\t$exon_intron3\t$$variant[9]\t$$variant[10]\t-1\t$gene4\t$exon_intron4\t$$variant[0]\t$$variant[1]\t$id[1]\t$mpd[1]\t$srd[1]\t$hm[1]\t$rp[1]\n";
	}
	# inter
	elsif ($$variant[0] =~ /inss/)
	{
		my ($gene1, $exon_intron1, $up_gene1, $up_distance1, $down_gene1, $down_distance1) = &mut_gene($$variant[5], $$variant[6]);
		my ($gene2, $exon_intron2, $up_gene2, $up_distance2, $down_gene2, $down_distance2) = &mut_gene($$variant[9], $$variant[10]);
		my ($gene3, $exon_intron3, $up_gene3, $up_distance3, $down_gene3, $down_distance3) = &mut_gene($$variant[5], $$variant[7]);
		my ($gene4, $exon_intron4, $up_gene4, $up_distance4, $down_gene4, $down_distance4) = &mut_gene($$variant[9], $$variant[11]);
		my @id = split ("\/", $$variant[2]);
		my @mpd = split ("\/", $$variant[3]);
		my @srd = split ("\/", $$variant[4]);
		my @hm = split ("\/", $$variant[13]);
		if ($$variant[5] lt $$variant[9])
		{
			$toprint = "$$variant[5]\t$$variant[6]\t1\t$gene1\t$exon_intron1\t$$variant[9]\t$$variant[10]\t-1\t$gene2\t$exon_intron2\t$$variant[0]\t$$variant[1]\t$id[0]\t$mpd[0]\t$srd[0]\t$hm[0]\t$rp[0]\n";
			$toprint .= "$$variant[5]\t$$variant[7]\t-1\t$gene3\t$exon_intron3\t$$variant[9]\t$$variant[11]\t1\t$gene4\t$exon_intron4\t$$variant[0]\t$$variant[1]\t$id[1]\t$mpd[1]\t$srd[1]\t$hm[1]\t$rp[1]\n";
		}
		else
		{
			$toprint = "$$variant[5]\t$$variant[7]\t1\t$gene1\t$exon_intron1\t$$variant[9]\t$$variant[11]\t-1\t$gene2\t$exon_intron2\t$$variant[0]\t$$variant[1]\t$id[0]\t$mpd[0]\t$srd[0]\t$hm[0]\t$rp[0]\n";
			$toprint .= "$$variant[5]\t$$variant[6]\t-1\t$gene3\t$exon_intron3\t$$variant[9]\t$$variant[10]\t1\t$gene4\t$exon_intron4\t$$variant[0]\t$$variant[1]\t$id[1]\t$mpd[1]\t$srd[1]\t$hm[1]\t$rp[1]\n";
		}
	}
	elsif ($$variant[0] =~ /inso/)
	{
		my ($gene1, $exon_intron1, $up_gene1, $up_distance1, $down_gene1, $down_distance1) = &mut_gene($$variant[5], $$variant[6]);
		my ($gene2, $exon_intron2, $up_gene2, $up_distance2, $down_gene2, $down_distance2) = &mut_gene($$variant[9], $$variant[11]);
		my ($gene3, $exon_intron3, $up_gene3, $up_distance3, $down_gene3, $down_distance3) = &mut_gene($$variant[5], $$variant[7]);
		my ($gene4, $exon_intron4, $up_gene4, $up_distance4, $down_gene4, $down_distance4) = &mut_gene($$variant[9], $$variant[10]);
		my @id = split ("\/", $$variant[2]);
		my @mpd = split ("\/", $$variant[3]);
		my @srd = split ("\/", $$variant[4]);
		my @hm = split ("\/", $$variant[13]);
		$toprint = "$$variant[5]\t$$variant[6]\t1\t$gene1\t$exon_intron1\t$$variant[9]\t$$variant[11]\t1\t$gene2\t$exon_intron2\t$$variant[0]\t$$variant[1]\t$id[0]\t$mpd[0]\t$srd[0]\t$hm[0]\t$rp[0]\n";
		$toprint .= "$$variant[5]\t$$variant[7]\t-1\t$gene3\t$exon_intron3\t$$variant[9]\t$$variant[10]\t-1\t$gene4\t$exon_intron4\t$$variant[0]\t$$variant[1]\t$id[1]\t$mpd[1]\t$srd[1]\t$hm[1]\t$rp[1]\n";
	}
	elsif ($$variant[0] eq 'transl_inter')
	{
		my ($gene1, $exon_intron1, $up_gene1, $up_distance1, $down_gene1, $down_distance1) = &mut_gene($$variant[5], $$variant[6]);
		my ($gene2, $exon_intron2, $up_gene2, $up_distance2, $down_gene2, $down_distance2) = &mut_gene($$variant[8], $$variant[9]);
		$toprint = "$$variant[5]\t$$variant[6]\t$$variant[7]\t$gene1\t$exon_intron1\t$$variant[8]\t$$variant[9]\t$$variant[10]\t$gene2\t$exon_intron2\t$$variant[0]\t$$variant[1]\t$$variant[2]\t$$variant[3]\t$$variant[4]\t$$variant[11]\t$rp[0]\n";
	}
		
	print OUT "$toprint";
}
close OUT;

my @variant;
open FILE, "<$mutgenefile";
while ($newline = <FILE>)
{
	chomp $newline;
	my @data = split (/\t/, $newline);
	push @variant, \@data;
}
close FILE;

open OUT, ">$fusiontypefile";
foreach (@variant)
{
	my $variant = $_;
	my $rp;
	if ($$variant[-1] =~ /_/)
	{
		$rp = pop(@$variant);
	}
	my $toprint;
	my ($type1, $type2, $type3);
	
	my ($partner1, $partner2);
	if ($$variant[3])
	{
		if ($$variant[2] == 1)
		{
			if ($exon{$$variant[3]}[0] eq '+')
			{
				$partner1 = '5';
			}
			else
			{
				$partner1 = '3';
			}
		}
		else
		{
			if ($exon{$$variant[3]}[0] eq '+')
			{
				$partner1 = '3';
			}
			else
			{
				$partner1 = '5';
			}
		}
	}
	if ($$variant[8])
	{
		if ($$variant[7] == 1)
		{
			if ($exon{$$variant[8]}[0] eq '+')
			{
				$partner2 = '5';
			}
			else
			{
				$partner2 = '3';
			}
		}
		else
		{
			if ($exon{$$variant[8]}[0] eq '+')
			{
				$partner2 = '3';
			}
			else
			{
				$partner2 = '5';
			}
		}
	}
	my $partner = $partner1.'-'.$partner2;
	# $type1
	if ($$variant[3])
	{
		if ($$variant[8])
		{
			$type1 = 'gene-gene';
		}
		else
		{
			$type1 = 'gene-intergenic';
		}
	}
	else
	{
		if ($$variant[8])
		{
			$type1 = 'gene-intergenic';
		}
		else
		{
			$type1 = 'intergenic-intergenic';
		}
	}
	
	# $type2 gene-gene
	if ($type1 eq 'gene-gene')
	{
		if ($$variant[2] == 1 and $$variant[7] == -1)
		{
			if ($exon2{$$variant[3]}[0][0] eq '+' and $exon2{$$variant[8]}[0][0] eq '-')
			{
				$type2 = 'head-head';
			}
			if ($exon2{$$variant[3]}[0][0] eq '+' and $exon2{$$variant[8]}[0][0] eq '+')
			{
				$type2 = 'head-tail';
			}
			if ($exon2{$$variant[3]}[0][0] eq '-' and $exon2{$$variant[8]}[0][0] eq '-')
			{
				$type2 = 'head-tail';
			}
			if ($exon2{$$variant[3]}[0][0] eq '-' and $exon2{$$variant[8]}[0][0] eq '+')
			{
				$type2 = 'tail-tail';
			}
		}
		if ($$variant[2] == 1 and $$variant[7] == 1)
		{
			if ($exon2{$$variant[3]}[0][0] eq '+' and $exon2{$$variant[8]}[0][0] eq '-')
			{
				$type2 = 'head-tail';
			}
			if ($exon2{$$variant[3]}[0][0] eq '+' and $exon2{$$variant[8]}[0][0] eq '+')
			{
				$type2 = 'head-head';
			}
			if ($exon2{$$variant[3]}[0][0] eq '-' and $exon2{$$variant[8]}[0][0] eq '-')
			{
				$type2 = 'tail-tail';
			}
			if ($exon2{$$variant[3]}[0][0] eq '-' and $exon2{$$variant[8]}[0][0] eq '+')
			{
				$type2 = 'head-tail';
			}
		}
		if ($$variant[2] == -1 and $$variant[7] == -1)
		{
			if ($exon2{$$variant[3]}[0][0] eq '+' and $exon2{$$variant[8]}[0][0] eq '-')
			{
				$type2 = 'head-tail';
			}
			if ($exon2{$$variant[3]}[0][0] eq '+' and $exon2{$$variant[8]}[0][0] eq '+')
			{
				$type2 = 'tail-tail';
			}
			if ($exon2{$$variant[3]}[0][0] eq '-' and $exon2{$$variant[8]}[0][0] eq '-')
			{
				$type2 = 'head-head';
			}
			if ($exon2{$$variant[3]}[0][0] eq '-' and $exon2{$$variant[8]}[0][0] eq '+')
			{
				$type2 = 'head-tail';
			}
		}
		if ($$variant[2] == -1 and $$variant[7] == 1)
		{
			if ($exon2{$$variant[3]}[0][0] eq '+' and $exon2{$$variant[8]}[0][0] eq '-')
			{
				$type2 = 'tail-tail';
			}
			if ($exon2{$$variant[3]}[0][0] eq '+' and $exon2{$$variant[8]}[0][0] eq '+')
			{
				$type2 = 'head-tail';
			}
			if ($exon2{$$variant[3]}[0][0] eq '-' and $exon2{$$variant[8]}[0][0] eq '-')
			{
				$type2 = 'head-tail';
			}
			if ($exon2{$$variant[3]}[0][0] eq '-' and $exon2{$$variant[8]}[0][0] eq '+')
			{
				$type2 = 'head-head';
			}
		}
	}
	
	# $type3 head-tail
	if ($type2 eq 'head-tail')
	{
		my ($ref_type, $ref_ei1, $ref_ei2) = &headtail($$variant[3], $$variant[1], $$variant[2], $$variant[8], $$variant[6], $$variant[7]);
		
		my $i = 0;
		$type3 = $$ref_type[$i];
		$$variant[4] = $$ref_ei1[$i];
		$$variant[9] = $$ref_ei2[$i];
		foreach (@{$ref_type})
		{
			unless ($_ eq 'unknown')
			{
				if ($type3 eq 'unknown')
				{
					$type3 = $_;
					$$variant[4] = $$ref_ei1[$i];
					$$variant[9] = $$ref_ei2[$i];
					last;
				}
			}
			$i++;
		}
		my $i = 0;
		foreach (@{$ref_type})
		{
			if ($_ eq 'in_frame')
			{
				$type3 = $_;
				$$variant[4] = $$ref_ei1[$i];
				$$variant[9] = $$ref_ei2[$i];
				last;
			}
			$i++;
		}
		my $i = 0;
		foreach (@{$ref_type})
		{
			if ($_ eq 'no_impact')
			{
				$type3 = $_;
				$$variant[4] = $$ref_ei1[$i];
				$$variant[9] = $$ref_ei2[$i];
				last;
			}
			$i++;
		}
		
	}
	$toprint = "$type1\t$type2\t$type3\t";
	$toprint .= join ("\t", @$variant);
	$toprint .= "\t".$partner;
	$toprint .= "\t$rp" if ($rp);
	print OUT "$toprint\n";
}
close OUT;
system "rm $mutgenefile";


sub mut_gene
{
	my $chr = shift;
	my $pos = shift;
	my ($gene, $exon_intron, $up_gene, $up_distance, $down_gene, $down_distance);
	
	foreach (@{$gene{$chr}})
	{
		if ($pos < $$_[1])
		{
			$down_distance = $$_[1] - $pos;
			$down_gene = $$_[0];
			$down_distance = int($down_distance/1000+0.5);
			$down_distance = '+'.$down_distance.'k';
			last;
		}
		elsif ($pos > $$_[2])
		{
			$up_distance = $pos - $$_[2];
			$up_gene = $$_[0];
			next;
		}
		elsif (&covered($pos, $pos, $$_[1], $$_[2]))
		{
			$gene = $$_[0];
			my $ei = 0; #exon id
			my $ii = 0; #intron id
			my $last_exon;
			if ($exon{$gene}[0] eq '+')
			{
				foreach (@{$exon{$gene}[1]})
				{
					my $current_exon = $_;
					$ei++;
					#print "$gene\t$$current_exon[0], $$current_exon[1]\n";
					if (&covered($pos, $pos, $$current_exon[0], $$current_exon[1]))
					{
						$exon_intron = 'E'.$ei;
						last;
					}
					else
					{
						if ($last_exon)
						{
							$ii++;
							if (&covered($pos, $pos, $$last_exon[1], $$current_exon[0]))
							{
								$exon_intron = 'I'.$ii;
								last;
							}
						}
					}
					$last_exon = $current_exon;
				}
			}
			else
			{
				my @reverse_exon = reverse (@{$exon{$gene}[1]});
				foreach (@reverse_exon)
				{
					my $current_exon = $_;
					$ei++;
					#print "$gene\t$$current_exon[0], $$current_exon[1]\n";
					if (&covered($pos, $pos, $$current_exon[0], $$current_exon[1]))
					{
						$exon_intron = 'E'.$ei;
						last;
					}
					else
					{
						if ($last_exon)
						{
							$ii++;
							if (&covered($pos, $pos, $$current_exon[1], $$last_exon[0]))
							{
								$exon_intron = 'I'.$ii;
								last;
							}
						}
					}
					$last_exon = $current_exon;
				}
			}
		}
	}
	$up_distance = int($up_distance/1000+0.5);
	$up_distance = '-'.$up_distance.'k';
	#print "$gene, $exon_intron, $up_gene, $up_distance, $down_gene, $down_distance\n";
	return ($gene, $exon_intron, $up_gene, $up_distance, $down_gene, $down_distance);
}

sub headtail
{
	my $gene1 = shift;
	my $pos1 = shift;
	my $ori1 = shift;
	my $gene2 = shift;
	my $pos2 = shift;
	my $ori2 = shift;
	
	my (@type, @exon_introna, @exon_intronb);
	
	foreach (@{$exon2{$gene1}})
	{
		my $genea = $_;
		my $exon_introna;
		
		my $eia = 0; #exon id
		my $iia = 0; #intron id
		my $last_exona;
		if ($$genea[0] eq '+')
		{
			foreach (@{$$genea[1]})
			{
				my $current_exona = $_;
				$eia++;
				#print "$gene\t$$current_exon[0], $$current_exon[1]\n";
				if (&covered($pos1, $pos1, $$current_exona[0], $$current_exona[1]))
				{
					$exon_introna = 'E'.$eia;
					last;
				}
				else
				{
					if ($last_exona)
					{
						$iia++;
						if (&covered($pos1, $pos1, $$last_exona[1], $$current_exona[0]))
						{
							$exon_introna = 'I'.$iia;
							last;
						}
					}
				}
				$last_exona = $current_exona;
			}
		}
		else
		{
			my @reverse_exona = reverse (@{$$genea[1]});
			foreach (@reverse_exona)
			{
				my $current_exona = $_;
				$eia++;
				#print "$gene\t$$current_exon[0], $$current_exon[1]\n";
				if (&covered($pos1, $pos1, $$current_exona[0], $$current_exona[1]))
				{
					$exon_introna = 'E'.$eia;
					last;
				}
				else
				{
					if ($last_exona)
					{
						$iia++;
						if (&covered($pos1, $pos1, $$current_exona[1], $$last_exona[0]))
						{
							$exon_introna = 'I'.$iia;
							last;
						}
					}
				}
				$last_exona = $current_exona;
			}
		}
		
		my $framea;
		if ($iia)
		{
			if ($$genea[0] eq '+')
			{
				$framea = $$genea[1][$iia][2];
			}
			else
			{
				$framea = $$genea[1][-$iia][2];
			}
		}
			
		
		foreach (@{$exon2{$gene2}})
		{
			my $geneb = $_;
			my $exon_intronb;
			
			my $eib = 0; #exon id
			my $iib = 0; #intron id
			my $last_exonb;
			if ($$geneb[0] eq '+')
			{
				foreach (@{$$geneb[1]})
				{
					my $current_exonb = $_;
					$eib++;
					#print "$gene\t$$current_exon[0], $$current_exon[1]\n";
					if (&covered($pos2, $pos2, $$current_exonb[0], $$current_exonb[1]))
					{
						$exon_intronb = 'E'.$eib;
						last;
					}
					else
					{
						if ($last_exonb)
						{
							$iib++;
							if (&covered($pos2, $pos2, $$last_exonb[1], $$current_exonb[0]))
							{
								$exon_intronb = 'I'.$iib;
								last;
							}
						}
					}
					$last_exonb = $current_exonb;
				}
			}
			else
			{
				my @reverse_exonb = reverse (@{$$geneb[1]});
				foreach (@reverse_exonb)
				{
					my $current_exonb = $_;
					$eib++;
					#print "$gene\t$$current_exon[0], $$current_exon[1]\n";
					if (&covered($pos2, $pos2, $$current_exonb[0], $$current_exonb[1]))
					{
						$exon_intronb = 'E'.$eib;
						last;
					}
					else
					{
						if ($last_exonb)
						{
							$iib++;
							if (&covered($pos2, $pos2, $$current_exonb[1], $$last_exonb[0]))
							{
								$exon_intronb = 'I'.$iib;
								last;
							}
						}
					}
					$last_exonb = $current_exonb;
				}
			}
			
			my $frameb;
			if (defined($iib))
			{
				if ($$geneb[0] eq '+')
				{
					$frameb = $$geneb[1][$iib][2];
				}
				else
				{
					$frameb = $$geneb[1][-$iib][2];
				}
			}
			
			my $type;
			next unless ($exon_introna and $exon_intronb);
			if ($exon_introna and $exon_intronb)
			{
				if ($framea == $frameb)
				{
					$type = 'in_frame';
				}
				else
				{
					$type = 'out_of_frame';
				}
			}
			if ($exon_introna =~ /E/ or $exon_intronb =~ /E/)
			{
				$type = 'unknown';
			}
			if (&covered($pos1, $pos1, $$genea[2], $$genea[4]) or &covered($pos1, $pos1, $$genea[3], $$genea[5]) or &covered($pos2, $pos2, $$geneb[2], $$geneb[4]) or &covered($pos2, $pos2, $$geneb[3], $$geneb[5]))
			{
				$type = 'UTR_swap';
			}
			if ($gene1 eq $gene2 and $ori1 ne $ori2 and $genea eq $geneb and $exon_introna and $exon_intronb and $exon_introna =~ /I/ and $iia == $iib)
			{
				$type = 'no_impact';
			}
			push @type, $type;
			push @exon_introna, $exon_introna;
			push @exon_intronb, $exon_intronb;
			#print "$gene1\t$exon_introna\t$framea\t$gene2\t$exon_intronb\t$frameb\t$type\n";
		}
	}
	
	my @out = (\@type, \@exon_introna, \@exon_intronb);
	return @out;
}

sub covered
{
	my $a1 = shift;
	my $a2 = shift;
	my $b1 = shift;
	my $b2 = shift;
	return 0 if ($a1 =~ /\D/ or $a2 =~ /\D/ or $b1 =~ /\D/ or $b2 =~ /\D/);
	($a1, $a2) = sort{ $a <=> $b } ($a1, $a2);
	($b1, $b2) = sort { $a <=> $b } ($b1, $b2);
	if (($a1 >= $b1 and $a1 <= $b2) or ($a2 >= $b1 and $a2 <= $b2) or ($b1 >= $a1 and $b1 <= $a2) or ($b2 >= $a1 and $b2 <= $a2))
	{
		return 1;
	}
}

sub print_usage
{
	die "Usage:\nperl ./scripts/fusions.pl [options]
	-i FILE	input file name, required
	-G FILE	reference gene annotation file name, required
	-h help\n";
}
//============================================================================
// Name        : ctea.cpp
// Author      : Kyu Park ( park@kyu.co )
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

//#include <iostream>
//#include <bwtaln.h>
#include "ReadClipper.hpp"

using namespace TeaQ;
using namespace std;

int main(int argc, char *argv[]) {
	cout << "ctea-170613-9-teaify0.4.awk\n";
	string program_name;
	if (argc > 1) {
		program_name = argv[1];
	}

	if (argc != 4) {
		cerr << "usage:\t" << argv[0] << " <bam_filename> <ref_fa_filename>\n";
		cerr << "\t ref family file (example: /home/el114/kyu/vincent/Tea/lib/assembly/human_youngTE_revisedPolyA.fa)\n";
		return 1;
	}
	string bam_filename = argv[2]; // string bam_file_name = "/home/el114/kyu/data/test.bam";
	string ref_fa_filename = argv[3]; //"/home/el114/kyu/data/bwa_idx/human_youngTE_revisedPolyA.fa";
	int64_t minimum_read_length = 5;
	int64_t minimum_base_gap = 2;
	string cap3_options = " -i 21 -j 31 -o 16 -s 251 -p 70 ";

	string f_contigs_file_name = bam_filename + "-ctea/04-contigs-f";
	string r_contigs_file_name = bam_filename + "-ctea/04-contigs-r";

	if ("ctea" == program_name) {
		ReadClipper::clip_filter_contiggen(
				bam_filename, ref_fa_filename,
				minimum_read_length, minimum_base_gap, cap3_options);
	} else if ("cteax" == program_name) {
		ReadClipper::clip_filter_pickLongestAsContig(
				bam_filename, ref_fa_filename,
				minimum_read_length, minimum_base_gap);
	} else if ("cteax03" == program_name) {
		ReadClipper::pickLongestAsContig(bam_filename, 'f', ref_fa_filename);
		ReadClipper::pickLongestAsContig(bam_filename, 'r', ref_fa_filename);

		ReadClipper::run_bwa(f_contigs_file_name, ref_fa_filename);
		ReadClipper::run_bwa(r_contigs_file_name, ref_fa_filename);

		ReadClipper::filter_family(bam_filename);

		ReadClipper::give_refName(bam_filename);

	} else if ("ctea04" == program_name || "cteax04" == program_name) {
		ReadClipper::run_bwa(f_contigs_file_name, ref_fa_filename);
		ReadClipper::run_bwa(r_contigs_file_name, ref_fa_filename);

		ReadClipper::filter_family(bam_filename);

		ReadClipper::give_refName(bam_filename);

	} else if ("ctea05" == program_name || "cteax05" == program_name) {
		ReadClipper::give_refName(bam_filename);
	}

	return 0;
}

//============================================================================
// Name        : TeaQ.cpp
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
	cout << "cTea-170531-3-bad-lexical-cast\n";
	string program_name;
	if (argc > 1) {
		program_name = argv[1];
	}

	if (argc != 4) {
		cerr << "usage:\t" << argv[0] << " <bam_file_name> <ref_file_name>\n";
		cerr << "\t ref family file (example: /home/el114/kyu/vincent/Tea/lib/assembly/human_youngTE_revisedPolyA.fa)\n";
		return 1;
	}
	string bam_file_name = argv[2]; // string bam_file_name = "/home/el114/kyu/data/test.bam";
	string ref_file_name = argv[3]; //"/home/el114/kyu/data/bwa_idx/human_youngTE_revisedPolyA.fa";
	int64_t minimum_read_length = 5;
	int64_t minimum_base_gap = 2;
	string cap3_options = "-i 21 -j 31 -o 16 -s 251 -p 70";

	if ("ctea" == program_name) {
		ReadClipper::clip_filter_contiggen(
				bam_file_name, ref_file_name,
				minimum_read_length, minimum_base_gap, cap3_options);
	} else if ("ctea-express" == program_name) {
		ReadClipper::clip_filter_pickLongestAsContig(
				bam_file_name, ref_file_name,
				minimum_read_length, minimum_base_gap);
	} else if ("ctea05" == program_name) {
		ReadClipper::give_refName(bam_file_name);
	} else if ("ctea-express05" == program_name) {
		ReadClipper::give_refName(bam_file_name);
	}

	return 0;
}

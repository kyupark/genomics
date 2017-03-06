//============================================================================
// Name        : TeaQ.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

//#include <iostream>
//#include <bwtaln.h>
#include "ReadClipper.hpp"
#include <cstdlib>

using namespace TeaQ;
using namespace std;

int main(int argc, char *argv[]) {
	if (argc != 2) {
		cerr << "usage:\t" << argv[0] << " <bam_file_name> <ref_file_name>" << endl;
		cerr << "\t ref family file (default: /home/el114/kyu/data/bwa_idx/human_youngTE_revisedPolyA.fa)" << endl;
		return 1;
	}
	string bam_file_name = argv[1]; // string bam_file_name = "/home/el114/kyu/data/test.bam";
	string ref_file_name = argv[2]; //"/home/el114/kyu/data/bwa_idx/human_youngTE_revisedPolyA.fa";
	int64_t minimum_read_length = 5;
	int64_t minimum_base_gap = 2;
	string cap3_options = "-i 21 -j 31 -o 16 -s 251 -p 70";

	ReadClipper::clip_filter_contiggen(
			bam_file_name,
			ref_file_name,
			minimum_read_length,
			minimum_base_gap,
			cap3_options);

	return 0;
}

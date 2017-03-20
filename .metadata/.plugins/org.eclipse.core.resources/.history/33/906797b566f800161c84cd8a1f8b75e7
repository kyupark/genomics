//============================================================================
// Name        : TeaQ.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

// #include <iostream>
//#include <bwtaln.h>
#include "ReadClipper.hpp"


using namespace TeaQ;
using namespace std;

int main() {
        string bam_file_name = "/home/el114/kyu/data/test.bam";
        //ReadClipper::clip_filter_sort(bam_file_name, 5, 2);

        ReadClipper::clip_filter_contiggen(
        		bam_file_name,
				"/home/el114/kyu/data/bwa_idx/human_youngTE_revisedPolyA.fa",
                5, // minimum_read_length = 5,
                2, // minimum_base_gap = 2,
                "-i 21 -j 31 -o 16 -s 251 -p 70");
        			// cap3_options = "-i 21 -j 31 -o 16 -s 251 -p 70");



        return 0;
}

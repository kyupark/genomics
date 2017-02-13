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
#include <bwtaln.h>


using namespace TeaQ;
using namespace std;

int main() {
        string bam_file_name = "/home/el114/kyu/data/test.bam";
        //ReadClipper::clip_filter_sort(bam_file_name, 5, 2);

        ReadClipper::clip_filter_contiggen(
        		bam_file_name,
                5,
                2,
                "-i 21 -j 31 -o 16 -s 251 -p 70");


        return 0;
}

//============================================================================
// Name        : TeaQ.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

// #include <iostream>
#include "ReadClipper.hpp"

using namespace TeaQ;
using namespace std;

int main() {
	string file_name = "/home/el114/kyu/data/test.bam";
//	ReadClipper::clip_reads(file_name, 5);
//	ReadClipper::filter_reads(file_name, 2);
	ReadClipper::generate_contigs(file_name);
	return 0;
}

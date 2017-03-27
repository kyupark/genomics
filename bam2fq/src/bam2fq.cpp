//============================================================================
// Name        : bam2fq.cpp
// Author      : Q
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "castle/OptionParser.hpp"
#include "castle/TimeChecker.hpp"
#include "sub_modules.hpp"

using namespace std;

int main(int argc, char **argv) {
	setvbuf(stdout, NULL, _IONBF, 0); // flush buffer immediately
	cat::Bam2Fastq b2fq;
	castle::OptionParser options(argc, argv);
	castle::TimeChecker checker;
	checker.setTarget("Bam2Fastq.main");
	checker.start();

	b2fq.set_option_parser(options);
	b2fq.process();

	cout << checker;

	return 0;
}

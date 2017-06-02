/*
 * ReadClipper.hpp
 *
 *  Created on: Feb 7, 2017
 *      Author: el114
 */

#ifndef READCLIPPER_HPP_
#define READCLIPPER_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <thread>
#include <algorithm>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <boost/lexical_cast.hpp>

namespace TeaQ {

using namespace std;
using namespace BamTools;

class ReadClipper {
public:
        ReadClipper();
        ~ReadClipper();

		set<string> program_names;

		static void give_refName(string& bam_file_name);

        static bool clip_filter_contiggen(
        	string& bam_file_name,
			string& ref_fa_file_name,
        	int64_t minimum_read_length = 5,
        	int64_t minimum_base_gap = 2,
        	string cap3_options = "-i 21 -j 31 -o 16 -s 251 -p 70");

        static bool clip_filter_pickLongestAsContig(
			string& bam_file_name,
			string& ref_fa_file_name,
			int64_t minimum_read_length = 5,
			int64_t minimum_base_gap = 2);

};


} /* namespace TeaQ */

#endif /* READCLIPPER_HPP_ */

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
#include <iterator>
#include <string>
#include <vector>
#include <limits>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <boost/lexical_cast.hpp>


namespace TeaQ {

using namespace std;
using namespace BamTools;

const bool FORWARD = true;
const bool REVERSE = false;

class ReadClipper {
public:
	ReadClipper();
	~ReadClipper();

	static bool clip_reads(string& file_name, int64_t minimum_read_length);
	static bool filter_reads(string& file_name, int64_t minimum_base_gap);

};


} /* namespace QT */

#endif /* READCLIPPER_HPP_ */

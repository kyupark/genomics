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
#include <thread>
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

        static bool clip_filter_contiggen(
                        string& bam_file_name,
                        int64_t minimum_read_length = 5,
                        int64_t minimum_base_gap = 2,
                        string cap3_options = "-i 21 -j 31 -o 16 -s 251 -p 70");

};


} /* namespace QT */

#endif /* READCLIPPER_HPP_ */

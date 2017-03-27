/*
 * OptionParser.hpp
 *
 *  Created on: Jun 4, 2015
 *      Author: Euncheon Lim @ Max-Planck-Institute for Developmental Biology
 *            : under the guidance of Prof. Detlef Weigel and Prof. Daniel Huson
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#ifndef OPTIONPARSER_HPP_
#define OPTIONPARSER_HPP_
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/cstdint.hpp>
#include <boost/format.hpp>
#include "TimeChecker.hpp"
#include "StringUtils.hpp"
#include "../third/gzstream.h"

namespace castle {
	using namespace std;
	class OptionParser {
		public:
			string input_filename;
			string output_path;

			int32_t n_cores;
		public:
			void show_help();
			OptionParser();
			OptionParser(int argc, char **argv);
			~OptionParser();
		public:
			bool is_complete();
			void expand_path(string& a_path);
			void expand_home_path(string& a_path);
			void print_options();
			string get_working_path(const string& a_path);
			string get_working_folder_path() const;

	};
} /* namespace castle */
#endif /* OPTIONPARSER_HPP_ */

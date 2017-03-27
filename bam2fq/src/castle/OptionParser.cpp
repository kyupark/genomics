/*
 * OptionParser.cpp
 *
 *  Created on: Jun 4, 2015
 *      Author: Euncheon Lim @ Max-Planck-Institute for Developmental Biology
 *            : under the guidance of Prof. Detlef Weigel and Prof. Daniel Huson
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#include "OptionParser.hpp"

namespace castle {

void OptionParser::show_help() {
//	cout << "Meerkat Preprocess: Parallelized bamreader (Ver. 0.0.1) \n"
//			"==============================================================================\n";
//	cout << "Syntax:    pmeerkat preprocess [options] <BAM input file>\n\n" << "Options:\n"
//			"  -b cutoff  Generate a blacklist of genomic positions with coverage\n"
//			"             greater than 'cutoff'.  The input BAM must be sorted if\n";

//	cout << "Contact:   Euncheon Lim <abysslover@gmail.com>\n";
	exit(1);
}

OptionParser::OptionParser() {
	TimeChecker checker;
	n_cores = checker.get_number_of_cores();
}

OptionParser::OptionParser(int argc, char **argv) {
	TimeChecker checker;
	n_cores = checker.get_number_of_cores();

	for (int i = 1; i < argc; ++i) {
		try {
			string argument = argv[i];
			if ("-i" == argument) {
				input_filename = string(argv[i + 1]);
			} else if ("-o" == argument) {
				output_path = string(argv[i + 1]);
			} else if ("-h" == argument) {
				show_help();
			}
		} catch (exception& ex) {
			cout << ex.what() << "\n";
			cout << argv[i] << "/" << argv[i + 1] << "\n";
			exit(1);
		}
	}

	if (!is_complete()) {
		show_help();
	}

//	if ("preprocess" == program_name) {
//		print_options();
//	} else {
//		read_isinfo();
//		read_point_black_lists();
//	}

}

void OptionParser::print_options() {
//	cout << "Options:\n   input BAM: " << fname << "\n   prefix: " << prefix << "\n";
//	if (coverage_cutoff > 0) {
//		cout << "   generating blacklist: coverage >= " << coverage_cutoff << " reads\n";
//	} else {
//		cout << "   no blacklist will be generated\n";
//	}
//	cout << "   max insert size: " << max_isize << " bases\n" << "   soft clip threshold: " << big_s_bps << " bases\n" << "   fragment size: " << frag_size << " bases\n" << "   insert size samples: ";
//	if (isize_samples < 0) {
//		cout << "all";
//	} else {
//		cout << isize_samples;
//	}
//	cout << " reads\n";
//	cout << "   blacklisted read groups: ";
//
//	for (auto it = rg_blacklist.begin(); it != rg_blacklist.end(); ++it) {
//		cout << *it << " ";
//	}
//	cout << "\n";
//	cout << "   trim reads: " << q << "\n   max N bps per read: " << n_cutoff << "\n";
//	cout << "   output files: \n      " << umfname << "\n      " << clipname << "\n      " << isinfoname << "\n      " << umrdistname << "\n      " << scrdistname << "\n      " << split1name << "\n      " << split2name << "\n";
//	if (coverage_cutoff > 0) {
//		cout << "      " << blistname << "\n";
//	}
//	cout << "      " << prefix << "/ (directory)" << "\n";
//
//	if (!no_prompt) {
//		cout << "Continue? (y/N) ";
//		char c;
//		cin >> noskipws >> c;
//		if (tolower(c) != 'y') {
//			cout << "Terminating...\n";
//			exit(1);
//		}
//	} else {
//		cout << "Forcing immediate execution.\n";
//	}
}

OptionParser::~OptionParser() {
}

bool OptionParser::is_complete() {
	if (!input_filename.empty() && !output_path.empty()) {
		return true;
	}
	return false;
}
void OptionParser::expand_home_path(string& a_path) {
	if (a_path.empty() || '~' != a_path[0]) {
		return;
	}
	char const* home = getenv("HOME");
	if (home || ((home = getenv("USERPROFILE")))) {
		a_path.replace(0, 1, home);
	} else {
		char const *hdrive = getenv("HOMEDRIVE"), *hpath = getenv("HOMEPATH");
		a_path.replace(0, 1, string(hdrive) + hpath);
	}
}

void OptionParser::expand_path(string& a_path) {
	boost::filesystem::path current_path(boost::filesystem::current_path());
	if ('/' != a_path[0] && '~' != a_path[0]) {
		a_path = current_path.string() + "/" + a_path;
	} else if ('~' == a_path[0]) {
		expand_home_path(a_path);
	}
}
string OptionParser::get_working_path(const string& a_path) {
	string copied_path(a_path);
	boost::filesystem::path current_path(boost::filesystem::current_path());
	string file_prefix;
	if ('/' != copied_path[0] && '~' != copied_path[0]) {
		file_prefix = current_path.string();
	} else if ('~' == copied_path[0]) {
		expand_home_path(copied_path);
	}
	file_prefix = copied_path.substr(0, copied_path.rfind('/'));
	return file_prefix;
}

} /* namespace castle */

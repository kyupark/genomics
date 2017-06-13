/*
 * filter-invisible-positives.cpp
 *
 *  Created on: Feb 24, 2017
 *      Author: el114
 */

#include "filter-invisible-positives.hpp"

using namespace FIP;
using namespace std;

int main(int argc, char *argv[]) {
	if (argc != 5) {
		cerr << "usage:\t" << argv[0] << " <wgs_file_format> <wes_file_format> <range from midpoint> <minimum_count>" << endl;
		return 1;
	}
	string wgs_file_name = argv[1];
	string wes_file_name = argv[2];
	int64_t range = 0;
	int64_t minimum_count = 0;

	istringstream ss3(argv[3]);
	if (!(ss3 >> range))
		cerr << "Invalid number " << argv[3] << '\n';

	istringstream ss4(argv[4]);
	if (!(ss4 >> minimum_count))
		cerr << "Invalid number " << argv[4] << '\n';



//  string wgs_file_name = "/home/el114/kyu/data/report-crcA56B_normal/WGS-10+read-set-intsct-Q-set--28uniq.WGS-10+read-set";
//	string wes_file_name = "/home/el174/vincent_lim/tea_test/exome/crcA56B_normal.ram.bam";

	filter_and_save(wgs_file_name, wes_file_name, range, minimum_count);


	return 0;
}

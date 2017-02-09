/*
 * ReadClipper.cpp
 *
 *  Created on: Feb 7, 2017
 *      Author: el114
 */

#include "ReadClipper.hpp"

namespace TeaQ {
ReadClipper::ReadClipper() {
	// TODO Auto-generated constructor stub
}
ReadClipper::~ReadClipper() {
	// TODO Auto-generated destructor stub
}

bool ReadClipper::clip_reads(string& file_name, int64_t minimum_read_length) {
	cout << "Clipping Reads" << endl;

	ifstream in(file_name, ifstream::binary);
	string line;

	BamReader reader;

	if ( !reader.Open(file_name)) {
		cerr << "Could not open input BAM file" << endl;
		return false;
	}

	BamAlignment al;
	string clipped_f_file_name = file_name + ".clipped_f";
	string clipped_r_file_name = file_name + ".clipped_r";
	ofstream clipped_f(clipped_f_file_name);
	cout << "Writing " << clipped_f_file_name << endl;
	ofstream clipped_r(clipped_r_file_name);
	cout << "Writing " << clipped_r_file_name << endl;

	while(reader.GetNextAlignment(al)) {
		if (al.IsDuplicate()) {
			continue;
		}

		int64_t refID = al.RefID;
		int64_t pos = al.Position;
		string bases = al.QueryBases;

//		int64_t bases_length = bases.length();
		stringstream ss_f;
		stringstream ss_r;

		int64_t cursor = 0;
		for(std::vector<int64_t>::size_type i = 0; i != al.CigarData.size(); ++i) {
			auto& a_cigar = al.CigarData[i];
			char& type = a_cigar.Type;
			int64_t len = a_cigar.Length;

//			stream << pos << "\t" << bases << endl;
			if (type == 'S' && len >= minimum_read_length) {
				// forward
				if (!al.IsReverseStrand()){
					if (i == 0){
						ss_f << refID << "\t" << pos + 1 << "\t" << bases.substr(cursor, len) << endl;
					}
				}
				// reverse
				else {
					if (i == al.CigarData.size()-1) {
						ss_r  << refID << "\t" << pos + 1 << "\t" << bases.substr(cursor, len) << endl;
					}
				}
			}
			else if (type == 'D') {
				continue;
			}

			cursor += len;
			pos += len;

			/*
			// checking because some reads do not have accurate cigar info
			// skip writing file
			if (cursor > bases_length) {
				break;
			}
			*/
		}

		//if (cursor == bases_length) {
			clipped_f << ss_f.str();
			clipped_r << ss_r.str();
		//}
	}

	cout << "Closing " << clipped_f_file_name << endl;
	clipped_f.close();

	cout << "Closing " << clipped_r_file_name << endl;
	clipped_r.close();

	cout << "Sorting " << clipped_f_file_name << endl;
	string sort_f_cmd = "sort -g -k1 -k2 " + clipped_f_file_name + " -o " + clipped_f_file_name; // + ".sorted";
	system(sort_f_cmd.c_str());

	cout << "Sorting " << clipped_r_file_name << endl;
	string sort_r_cmd = "sort -g -k1 -k2 " + clipped_r_file_name + " -o " + clipped_r_file_name; // + ".sorted";
	system(sort_r_cmd.c_str());


	return true;
}

bool ReadClipper::filter_reads(string& file_name, int64_t minimum_base_gap) {
	cout << "Filtering Reads" << endl;

	ifstream in_f(file_name + ".clipped_f", ifstream::binary);
	ifstream in_r(file_name + ".clipped_r", ifstream::binary);

	string line;

	string al_prev_refID;
	int64_t al_prev_pos;
	string al_prev_bases;

	int64_t al_current_refID_int;
	string al_current_refID;
	int64_t al_current_pos;
	string al_current_bases;

	string filtered_f_file_name = file_name + ".filtered_f";
	string filtered_r_file_name = file_name + ".filtered_r";
	ofstream filtered_f(filtered_f_file_name);
	cout << "Writing " << filtered_f_file_name << endl;
	ofstream filtered_r(filtered_r_file_name);
	cout << "Writing " << filtered_r_file_name << endl;

	bool is_group = false;

	while(getline(in_f, line)){

		stringstream linestream(line);
		string value;

		//stringstream ss;

		if(!line.empty()){

			getline(linestream, value, '\t');
			al_current_refID_int = boost::lexical_cast<int64_t>(value);
			++al_current_refID_int;

			if (al_current_refID_int == 23) {
				al_current_refID = "chrX";
			}
			else if (al_current_refID_int == 24) {
				al_current_refID = "chrY";
			}
			else {
				al_current_refID = "chr" + to_string(al_current_refID_int);
			}

			getline(linestream, value, '\t');
			al_current_pos = boost::lexical_cast<int64_t>(value);
			getline(linestream, value, '\t');
			al_current_bases = value;

			if (al_current_refID == al_prev_refID
					&& al_current_pos - al_prev_pos <= minimum_base_gap) {
				filtered_f << ">" << al_prev_refID << ":" << al_prev_pos << endl;
				filtered_f << al_prev_bases << endl;
				is_group = true;
			}
			else {
				if (is_group) {
					filtered_f << ">" << al_prev_refID << ":" << al_prev_pos << endl;
					filtered_f << al_prev_bases << endl;
					filtered_f << endl;
				}
				is_group = false;
			}

			al_prev_refID = al_current_refID;
			al_prev_pos = al_current_pos;
			al_prev_bases = al_current_bases;
		}
	}

	al_prev_refID = "";
	al_prev_pos = 0;
	al_prev_bases = "";

	al_current_refID = "";
	al_current_pos = 0;
	al_current_bases = "";

	while(getline(in_r, line)){

			stringstream linestream(line);
			string value;

			//stringstream ss;

			if(!line.empty()){

				getline(linestream, value, '\t');
				al_current_refID_int = boost::lexical_cast<int64_t>(value);
				++al_current_refID_int;

				if (al_current_refID_int == 23) {
					al_current_refID = "chrX";
				}
				else if (al_current_refID_int == 24) {
					al_current_refID = "chrY";
				}
				else {
					al_current_refID = "chr" + to_string(al_current_refID_int);
				}

				getline(linestream, value, '\t');
				al_current_pos = boost::lexical_cast<int64_t>(value);
				getline(linestream, value, '\t');
				al_current_bases = value;

				if (al_current_refID == al_prev_refID
						&& al_current_pos - al_prev_pos <= minimum_base_gap) {
					filtered_r << ">" << al_prev_refID << ":" << al_prev_pos << endl;
					filtered_r << al_prev_bases << endl;

					is_group = true;
				}
				else {
					if (is_group) {
						filtered_r << ">" << al_prev_refID << ":" << al_prev_pos << endl;
						filtered_r << al_prev_bases << endl;
						filtered_r << endl;
					}
					is_group = false;
				}

				al_prev_refID = al_current_refID;
				al_prev_pos = al_current_pos;
				al_prev_bases = al_current_bases;
			}
		}

	filtered_f.close();
	cout << "Closing " << filtered_f_file_name << endl;
	filtered_r.close();
	cout << "Closing " << filtered_r_file_name << endl;

	return true;
}


bool ReadClipper::generate_contigs(string& file_name, string options) {
	ifstream in_f(file_name + ".filtered_f", ifstream::binary);
	string line;
	stringstream ss;
	string cmd_cap3;
	string temp_fa_log_file_name;

	string cap3_output_dir = file_name + ".cap3_output";
	string cmd_rmdir = "rm -rf " + cap3_output_dir;
	system(cmd_rmdir.c_str());
	string cmd_mkdir = "mkdir " + cap3_output_dir;
	system(cmd_mkdir.c_str());

	string temp_fa_file_prefix = cap3_output_dir + "/temp.fa";
	string temp_fa_file_name = temp_fa_file_prefix + ".0";
	ofstream temp_fa(temp_fa_file_name);

	string contigs_file_name = file_name + ".contigs";
	ofstream contigs(contigs_file_name);

	int64_t count = 0;
	while(getline(in_f, line)){
		if (!line.empty()) {
			temp_fa << line << endl;
		}
		else {
			temp_fa_file_name = temp_fa_file_prefix + "." + to_string(count);
			temp_fa_log_file_name = temp_fa_file_name + ".log";

			temp_fa.close();
			cmd_cap3 = "cap3 " + temp_fa_file_name + " " + options +
					" > " + temp_fa_log_file_name;
			cout << cmd_cap3 << endl;
			system(cmd_cap3.c_str());

			ifstream in(temp_fa_log_file_name, ifstream::binary);
			string log_line;
			int64_t cursor = 0;
			int64_t no_of_reads = 0;
			if (getline(in, log_line)) {
				while(getline(in, log_line)) {
					if (log_line.find("***") == -1 && cursor == 1) {
						++no_of_reads;
					}
					if (log_line.find("***") == -1 && cursor < 2) {
						++cursor;
						continue;
					}
					else {
						getline(in, log_line);
						break;
					}
				}

				no_of_reads -= 2;

				int64_t colon;
				int64_t plus_sign;
				string ref;
				string pos;
				stringstream ss_contigs;

				for (int64_t i = 0; i != no_of_reads; ++i) {
					getline(in, log_line);
					colon = log_line.find(":");

					if (i == 0) {
						cout << log_line << endl << " colon:" << colon;
						ref = log_line.substr(0, colon);
						contigs << ref << ":";
					}
					plus_sign = log_line.find("+");
					cout << " +:" << plus_sign << endl;
					pos = log_line.substr(colon, plus_sign);

					if (i != 0) {
						contigs << ",";
					}
					contigs << pos;
				}
				break;
			}

			++count;
			temp_fa_file_name = cap3_output_dir + "/temp.fa" + "." + to_string(count);
			temp_fa.open(temp_fa_file_name);
		}
	}
	contigs.close();

	return true;
}

} /* namespace QT */



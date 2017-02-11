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

//bool ReadClipper::clip_sort_filter(string& bam_file_name,
//		int64_t minimum_read_length = 5,
//		int64_t minimum_base_gap = 2,
//		string cap3_options = "-i 21 -j 31 -o 16 -s 251 -p 70") {
//
//	clip_reads(string& file_name, int64_t minimum_read_length);
//	filter_reads(string& file_name, int64_t minimum_base_gap);
//	generate_contigs(string& file_name, string cap3_options);
//
//}

bool ReadClipper::clip_reads(string& bam_file_name, int64_t minimum_read_length) {
	cout << "Clipping Reads" << endl;

	BamReader reader;
	if ( !reader.Open(bam_file_name)) {
		cerr << "Could not open input BAM file" << endl;
		return false;
	}

	BamAlignment al;
	auto& refData = reader.GetReferenceData();
	string f_clipped_file_name = bam_file_name + ".f_clipped";
	string r_clipped_file_name = bam_file_name + ".r_clipped";
	ofstream f_clipped(f_clipped_file_name);
	cout << "Writing " << f_clipped_file_name << endl;
	ofstream r_clipped(r_clipped_file_name);
	cout << "Writing " << r_clipped_file_name << endl;

//	int count = 0;
	while(reader.GetNextAlignment(al)) {
//		++count;
//		if (count == 50) {
//			break;
//		}
		if (al.IsDuplicate()) {
			continue;
		}

		string refID = to_string(al.RefID);
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

		}

		cout << ss_f.str();

		f_clipped << ss_f.str();
		r_clipped << ss_r.str();

	}

	cout << "Closing " << f_clipped_file_name << endl;
	f_clipped.close();

	cout << "Closing " << r_clipped_file_name << endl;
	r_clipped.close();

	cout << "Sorting " << f_clipped_file_name << endl;
	string sort_f_cmd = "sort -g -k1 -k2 " + f_clipped_file_name + " -o " + f_clipped_file_name;
	system(sort_f_cmd.c_str());

	cout << "Sorting " << r_clipped_file_name << endl;
	string sort_r_cmd = "sort -g -k1 -k2 " + r_clipped_file_name + " -o " + r_clipped_file_name;
	system(sort_r_cmd.c_str());


	return true;
}

bool ReadClipper::filter_reads(string& bam_file_name, int64_t minimum_base_gap) {
	cout << "Filtering Reads" << endl;

	ifstream in_f(bam_file_name + ".f_clipped", ifstream::binary);
	ifstream in_r(bam_file_name + ".r_clipped", ifstream::binary);

	string line;

	string al_prev_refName;
	int64_t al_prev_pos;
	string al_prev_bases;

	string al_current_refName;
	int64_t al_current_pos;
	string al_current_bases;

	string filtered_f_file_name = bam_file_name + ".f_filtered";
	string filtered_r_file_name = bam_file_name + ".r_filtered";
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
			al_current_refName = value;
			getline(linestream, value, '\t');
			al_current_pos = boost::lexical_cast<int64_t>(value);
			getline(linestream, value, '\t');
			al_current_bases = value;

			if (al_current_refName == al_prev_refName
					&& al_current_pos - al_prev_pos <= minimum_base_gap) {
				filtered_f  << al_prev_refName << "\t" << al_prev_pos << "\t" << al_prev_bases << endl;
				is_group = true;
			}
			else {
				if (is_group) {
					filtered_f << al_prev_refName << "\t" << al_prev_pos << "\t" << al_prev_bases << endl << endl;
				}
				is_group = false;
			}

			al_prev_refName = al_current_refName;
			al_prev_pos = al_current_pos;
			al_prev_bases = al_current_bases;
		}
	}

	al_prev_refName = "";
	al_prev_pos = 0;
	al_prev_bases = "";

	al_current_refName = "";
	al_current_pos = 0;
	al_current_bases = "";

	filtered_f.close();
	cout << "Closing " << filtered_f_file_name << endl;
	filtered_r.close();
	cout << "Closing " << filtered_r_file_name << endl;

	return true;
}

bool ReadClipper::generate_contigs(string& bam_file_name, string cap3_options) {
	ifstream in_f(bam_file_name + ".f_filtered", ifstream::binary);
	string line;

	string cmd_cap3;
	string log_file_name;
	string cap_info_file_name;

	string cap3_output_dir = bam_file_name + ".cap3_output";
	string cmd_rmdir = "rm -rf " + cap3_output_dir;
	system(cmd_rmdir.c_str());
	string cmd_mkdir = "mkdir " + cap3_output_dir;
	system(cmd_mkdir.c_str());

	string temp_fa_file_prefix = cap3_output_dir + "/temp.fa";
	string temp_fa_file_name = temp_fa_file_prefix + ".0";
	ofstream temp_fa(temp_fa_file_name);

	string contigs_file_name = bam_file_name + ".f_contigs";
	ofstream contigs(contigs_file_name);

	string tea_file_name = bam_file_name + ".f_tea";
	ofstream tea(tea_file_name);

	int64_t count = 0;
	string chr;
	int64_t cpos = 0;
	int64_t pos;
	string pos_combined;
	string base;
	string base_combined;
	string contig;
	int64_t no_of_reads = 0;

	while(getline(in_f, line)){
		if (!line.empty()) {
			++no_of_reads;

			int64_t chr_index = line.find_first_of("\t");
			int64_t base_index = line.find_last_of("\t");

			chr = line.substr(0, chr_index);
			pos = boost::lexical_cast<int64_t>(line.substr(chr_index+1, base_index-chr_index-1));
			base = line.substr(base_index+1);
			temp_fa << ">" << chr << ":" << pos << endl;
			temp_fa << base << endl;

			if (cpos == 0){
				cpos = pos;
			}
			else if (cpos < pos) {
				cpos = pos;
			}

			if (pos_combined.length() == 0) {
				pos_combined = to_string(pos);
				base_combined = base;
			}
			else {
				pos_combined += "," + to_string(pos);
				base_combined += "," + base;
			}
		}
		else {
			//stringstream ss;

			temp_fa_file_name = temp_fa_file_prefix + "." + to_string(count);
			log_file_name = temp_fa_file_name + ".log";
			string cap_contigs_file_name = temp_fa_file_name + ".cap.contigs";
			string cap_singlets_file_name = temp_fa_file_name + ".cap.singlets";

			temp_fa.close();
			cmd_cap3 = "cap3 " + temp_fa_file_name + " " + cap3_options +
					" > " + log_file_name;
			cout << cmd_cap3 << endl;
			system(cmd_cap3.c_str());

			int64_t s_chr_index = 0;
			int64_t s_pos_index = 0;

			string line_contig;
			ifstream in_contigs(cap_contigs_file_name, ifstream::binary);
			getline(in_contigs, line_contig);

			// TODO need to handle multi line singlets
			// if contigs file is empty, check singlets file and get the longest
			if (line_contig.empty()) {
				string longest_singlet;
				string line_singlets;
				ifstream in_singlets(cap_singlets_file_name, ifstream::binary);
				if (!getline(in_singlets, line_singlets)) {
					break;
				}

				s_chr_index = line_singlets.find(">");
				s_pos_index = line_singlets.find(":");
				if (s_chr_index == string::npos
						|| s_pos_index == string::npos) {
					break;
				}

				chr = line_singlets.substr(s_chr_index+1, s_pos_index-s_chr_index-1);
				pos = boost::lexical_cast<int64_t>(line_singlets.substr(s_pos_index+1));

				while (getline(in_singlets, line_singlets)) {
					if (line_singlets.length() > longest_singlet.length()) {
						longest_singlet = line_singlets;
					}

					if (!getline(in_singlets, line_singlets)) {
						break;
					}

					s_chr_index = line_singlets.find(">");
					s_pos_index = line_singlets.find(":");
					if (s_chr_index == string::npos
							|| s_pos_index == string::npos) {
						break;
					}

					chr = line_singlets.substr(s_chr_index+1, s_pos_index-s_chr_index-1);
					pos = boost::lexical_cast<int64_t>(line_singlets.substr(s_pos_index+1));

				}
				contig = longest_singlet;
			}
			// if contigs file is not empty,
			else {
				getline(in_contigs, line_contig);
				contig = line_contig;

				// append multiline contig
				while(getline(in_contigs, line_contig)) {
					if (line_contig.find(">") == string::npos) {
						contig += line_contig;
					}
					// more than one contig
					else {
						cerr << "More than one contig: " << cap_contigs_file_name << endl;
						break;
					}
				}
			}

			tea << chr << "\t" << cpos << "\t" << no_of_reads << "\t" << contig << "\t";
			tea << pos_combined << "\t" << base_combined << endl;

			contigs << ">" << chr << ":" << cpos << endl;
			contigs << contig << endl;


			string rm_temp_file_cmd = "rm -rf " + temp_fa_file_name+ "*";
			system(rm_temp_file_cmd.c_str());

			++count;
			temp_fa_file_name = cap3_output_dir + "/temp.fa" + "." + to_string(count);
			temp_fa.open(temp_fa_file_name);

			// reinitializing variables
			chr = "";
			cpos = 0;
			pos = 0;
			pos_combined = "";
			base = "";
			base_combined = "";
			contig = "";
			no_of_reads = 0;
		}
	}
	contigs.close();

	return true;
}

//bwa aln bwa_idx/human_youngTE_revisedPolyA.fa test.bam.contigs > test.bam.contigs.sai
//bwa samse -n 100 bwa_idx/human_youngTE_revisedPolyA.fa test.bam.contigs.sai test.bam.contigs  > test.bam.contigs.bam

//bwa mem bwa_idx/human_youngTE_revisedPolyA.fa test.bam.contigs > test.bam.contigs.mem.bam



} /* namespace QT */



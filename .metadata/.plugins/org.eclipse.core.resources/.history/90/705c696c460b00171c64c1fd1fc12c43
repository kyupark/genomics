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

void sort(string clipped_file_name){
	string cmd_sort = "sort -g -k1 -k2 " + clipped_file_name + " -o " + clipped_file_name + "-sorted";
	system(cmd_sort.c_str());
	cout << "Sorted " << clipped_file_name << endl;
}

void rm_temp_all(string temp_file_name) {
	string cmd_rm_temp_all = "rm -rf " + temp_file_name+ "*";
	system(cmd_rm_temp_all.c_str());
}

void rm_temp(string temp_file_name) {
	string cmd_rm_temp = "rm -rf " + temp_file_name;
	system(cmd_rm_temp.c_str());
}

void mkdir(string bam_file_name) {
	string cmd_rm_debug_dir = "rm -rf " + bam_file_name + "-tea-debug";
	system(cmd_rm_debug_dir.c_str());
	string cmd_mk_debug_dir = "mkdir " + bam_file_name + "-tea-debug";
	system(cmd_mk_debug_dir.c_str());
}

void do_join(std::thread& t) {
    t.join();
}

void join_all(std::vector<std::thread>& tv) {
    for(auto& t : tv) {
    	do_join(t);
    }
}

void get_refName(string& bam_file_name) {
	BamReader reader;
	if ( !reader.Open(bam_file_name)) {
		cerr << "Could not open input BAM file" << endl;
		//return false;
	}

	auto& refData = reader.GetReferenceData();
	string refName_file_name = bam_file_name + "-tea-debug/00-refName";
	ofstream refName(refName_file_name);

	for (int64_t i=0; i < refData.size(); ++i) {
		refName << i << "\t" << refData[i].RefName << endl;
	}
}

void clip(string& bam_file_name, int64_t minimum_read_length) {
	cout << "Clipping Reads" << endl;

	BamReader reader;
	if ( !reader.Open(bam_file_name)) {
		cerr << "Could not open input BAM file" << endl;
		//return false;
	}

	BamAlignment al;
	string f_clipped_file_name = bam_file_name + "-tea-debug/01-clipped-f";
	string r_clipped_file_name = bam_file_name + "-tea-debug/01-clipped-r";
	ofstream f_clipped(f_clipped_file_name);
	cout << "Writing " << f_clipped_file_name << endl;
	ofstream r_clipped(r_clipped_file_name);
	cout << "Writing " << r_clipped_file_name << endl;

	while(reader.GetNextAlignment(al)) {
		if (al.IsDuplicate()) {
			continue;
		}

		string refID = to_string(al.RefID);
		int64_t pos = al.Position;
		string bases = al.QueryBases;

		stringstream ss_f;
		stringstream ss_r;

		int64_t cursor = 0;
		for(std::vector<int64_t>::size_type i = 0; i != al.CigarData.size(); ++i) {
			auto& a_cigar = al.CigarData[i];
			char& type = a_cigar.Type;
			int64_t len = a_cigar.Length;

			if (type == 'D' || type == 'd') {
				pos += len;
			}
			else if ( type == 'H' || type == 'h') {
				if (i == 0) {
					if (!al.IsReverseStrand() && len >= minimum_read_length) {
						ss_f << refID << "\t" << pos << "\t" << endl;
					}
					else {
						cursor += len;
						continue;
					}
				}
				else if (i == al.CigarData.size()-1) {
					if (al.IsReverseStrand() && len >= minimum_read_length) {
						ss_r << refID << "\t" << pos+1 << "\t" << endl;
					}
					else {
						cursor += len;
						continue;
					}
				}
				pos += len;
			}
			else if (type == 'S' | type == 's') {
				if (i == 0) {
					if (!al.IsReverseStrand() && len >= minimum_read_length) {
						ss_f << refID << "\t" << pos << "\t" << bases.substr(cursor, len) << endl;
					}
					else {
						cursor += len;
						continue;
					}
				}
				else if (i == al.CigarData.size()-1) {
					if (al.IsReverseStrand() && len >= minimum_read_length) {
						ss_r << refID << "\t" << pos+1 << "\t" << bases.substr(cursor, len) << endl;
					}
					else {
						cursor += len;
						continue;
					}
				}
				cursor += len;
				pos += len;
			}
			else {
				cursor += len;
				pos += len;
			}

		}

		f_clipped << ss_f.str();
		r_clipped << ss_r.str();

	}


	cout << "Closing " << f_clipped_file_name << endl;
	f_clipped.close();

	cout << "Closing " << r_clipped_file_name << endl;
	r_clipped.close();

	vector<std::thread> threads;

	threads.push_back(thread(sort, f_clipped_file_name));
	threads.push_back(thread(sort, r_clipped_file_name));

	join_all(threads);
}


void filter(string bam_file_name, char f_or_r, int64_t minimum_base_gap) {

	ifstream in(bam_file_name + "-tea-debug/01-clipped-" + f_or_r + "-sorted", ifstream::binary);

	string line;

	string al_prev_refID;
	int64_t al_prev_pos;
	string al_prev_bases;

	string al_current_refID;
	int64_t al_current_pos;
	string al_current_bases;

	string filtered_file_name = bam_file_name + "-tea-debug/02-filtered-" + f_or_r ;
	ofstream filtered(filtered_file_name);
	cout << "Writing " << filtered_file_name << endl;

	bool is_previous_space = true;
	int64_t count = 0;

	while(getline(in, line)){

		stringstream linestream(line);
		string value;

		if(!line.empty()){

			getline(linestream, value, '\t');
			al_current_refID = value;
			getline(linestream, value, '\t');
			al_current_pos = boost::lexical_cast<int64_t>(value);
			getline(linestream, value, '\t');
			al_current_bases = value;

			if (count == 0) {
				++count;
			}
			else if (count == 1) {
				if (al_current_refID == al_prev_refID
						&& abs(al_current_pos-al_prev_pos) <= minimum_base_gap) {
					filtered << al_prev_refID << "\t" << al_prev_pos << "\t" << al_prev_bases << endl;
					filtered << al_current_refID << "\t" << al_current_pos << "\t" << al_current_bases << endl;
					is_previous_space = false;
					++count;
				}
				else {
					count = 1;
					if (!is_previous_space) {
						filtered << endl;
						is_previous_space = true;
					}
				}
			}
			else if (count > 1) {
				if (al_current_refID == al_prev_refID
						&& abs(al_current_pos-al_prev_pos) <= minimum_base_gap) {
					filtered << al_current_refID << "\t" << al_current_pos << "\t" << al_current_bases << endl;
					is_previous_space = false;
					++count;
				}
				else {
					count = 1;
					if (!is_previous_space) {
						filtered << endl;
						is_previous_space = true;
					}
				}
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

	filtered.close();
	cout << "Closing " << filtered_file_name << endl;

}


void run_filter(string& bam_file_name, int64_t minimum_base_gap) {
	cout << "Filtering Reads" << endl;

	vector<std::thread> threads;

	threads.push_back(thread(filter, bam_file_name, 'f', minimum_base_gap));
	threads.push_back(thread(filter, bam_file_name, 'r', minimum_base_gap));

	join_all(threads);
}


void bwa_aln_samse(string& contigs_file_name, string& ref_fa_file_name) {
	//bwa aln bwa_idx/human_youngTE_revisedPolyA.fa test.bam.contigs > test.bam.contigs.sai
	//bwa samse -n 100 bwa_idx/human_youngTE_revisedPolyA.fa test.bam.contigs.sai test.bam.contigs  > test.bam.contigs.bam

	string cmd_load_bwa ="module load seq/bwa/0.6.2";
	system(cmd_load_bwa.c_str());

	string bwa_sai_file_name = contigs_file_name + ".sai";
	string cmd_bwa_aln =
			"bwa aln "
			+ ref_fa_file_name + " "
			+ contigs_file_name
			+ " > " + bwa_sai_file_name;
	system(cmd_bwa_aln.c_str());

//	vector<std::thread> threads;
//	threads.push_back(thread(rm_temp, contigs_file_name));

	string cmd_bwa_samse =
			"bwa samse -n 100 "
			+ ref_fa_file_name + " "
			+ bwa_sai_file_name + " "
			+ contigs_file_name
			+ " > " + contigs_file_name + ".samse-refid";
	system(cmd_bwa_samse.c_str());

//	threads.push_back(thread(rm_temp, bwa_sai_file_name));
//	join_all(threads);
}


void contiggen(string bam_file_name, char f_or_r, string cap3_options, string ref_fa_file_name) {
	vector<std::thread> threads;

	ifstream in(bam_file_name + "-tea-debug/02-filtered-" + f_or_r , ifstream::binary);
	string line;

	string cmd_cap3;
	string log_file_name;
	string cap_info_file_name;

	string cap3_output_dir = bam_file_name + "-tea-debug/03-cap3-output/";

	int64_t count = 0;
	string temp_fa_file_prefix = cap3_output_dir + "temp.fa." + f_or_r;
	string temp_fa_file_name = temp_fa_file_prefix + "." + to_string(count);
	ofstream temp_fa(temp_fa_file_name);

	string contigs_file_name = bam_file_name + "-tea-debug/04-contigs-" + f_or_r;
	ofstream contigs(contigs_file_name);

	string chr;
	int64_t refID;
	int64_t cpos = 0;
	int64_t pos;
	string pos_combined;
	string base;
	string base_combined;
	string contig;
	int64_t no_of_reads = 0;

	while(getline(in, line)){
		if (!line.empty()) {
			++no_of_reads;

			int64_t chr_index = line.find_first_of("\t");
			int64_t base_index = line.find_last_of("\t");

			chr = line.substr(0, chr_index);
			pos = boost::lexical_cast<int64_t>(line.substr(chr_index+1, base_index-chr_index-1));
			base = line.substr(base_index+1);

			if (base == "") continue;

			temp_fa << ">" << chr << ":" << pos << endl;
			temp_fa << base << endl;

			if (cpos == 0){
				cpos = pos;
			}
			else if (cpos > pos) {
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
			temp_fa_file_name = temp_fa_file_prefix + "." + to_string(count);
			log_file_name = temp_fa_file_name + ".log";
			string cap_contigs_file_name = temp_fa_file_name + ".cap.contigs";
			string cap_singlets_file_name = temp_fa_file_name + ".cap.singlets";

			temp_fa.close();

			// In case, there were only Hard clips(empty bases)
			if (base_combined != "") {

				cmd_cap3 = "cap3 " + temp_fa_file_name + " " + cap3_options +
						" > " + log_file_name;
				// cout << cmd_cap3 << endl;
				system(cmd_cap3.c_str());

				int64_t s_chr_index = 0;
				int64_t s_pos_index = 0;

				string line_contigs;
				ifstream in_contigs(cap_contigs_file_name, ifstream::binary);
				getline(in_contigs, line_contigs);

				// if contigs file is empty, check singlets file and get the longest
				if (line_contigs.empty()) {
					string longest_singlet="";
					string line_singlets="";
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
					getline(in_contigs, line_contigs);
					contig = line_contigs;

					// append multiline contig
					while(getline(in_contigs, line_contigs)) {
						if (line_contigs.find(">") == string::npos) {
							contig += line_contigs;
						}
						// TODO handle more than one contig
						else {
							cerr << "More than one contig: " << cap_contigs_file_name << endl;
							break;
						}
					}
				}

				contigs << ">" << chr << ";" << cpos << ";" << f_or_r << ";" << no_of_reads << ";" << pos_combined << ";" << base_combined << endl;
				contigs << contig << endl;

				rm_temp_all(temp_fa_file_name);

			}

			++count;
			temp_fa_file_name = temp_fa_file_prefix + "." + to_string(count);
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

//	threads.push_back(thread(rm_temp, bam_file_name+ ".clipped_" + f_or_r));
//	threads.push_back(thread(rm_temp, bam_file_name+ ".filtered_" + f_or_r));

	bwa_aln_samse(contigs_file_name, ref_fa_file_name);

//	join_all(threads);
}



void run_contiggen(string& bam_file_name, string& cap3_options, string& ref_fa_file_name) {
	cout << "Generating Reads" << endl;

	string cap3_output_dir = bam_file_name + "-tea-debug/03-cap3-output";
	string cmd_rmdir = "rm -rf " + cap3_output_dir;
	system(cmd_rmdir.c_str());
	string cmd_mkdir = "mkdir " + cap3_output_dir;
	system(cmd_mkdir.c_str());

	vector<std::thread> threads;

	threads.push_back(thread(contiggen, bam_file_name, 'f', cap3_options, ref_fa_file_name));
	threads.push_back(thread(contiggen, bam_file_name, 'r', cap3_options, ref_fa_file_name));

	join_all(threads);
}


void bwa_mem(string& contigs_file_name, string& ref_fa_file_name) {
	//bwa mem bwa_idx/human_youngTE_revisedPolyA.fa test.bam.contigs > test.bam.contigs.mem.bam

	string cmd_load_bwa ="module load seq/bwa/0.7.8";
	system(cmd_load_bwa.c_str());

	string cmd_bwa_mem =
			"bwa mem "
			+ ref_fa_file_name + " "
			+ contigs_file_name
			+ " > " + contigs_file_name + ".mem.bam";
	system(cmd_bwa_mem.c_str());
}

void run_bwa(string& contigs_file_name, string& ref_fa_file_name) {
	bwa_aln_samse(contigs_file_name, ref_fa_file_name);
	//bwa_mem(contigs_file_name, ref_fa_file_name);

}


void filter_family(string& bam_file_name) {
	string contigs_f_bam_file_name = bam_file_name + "-tea-debug/04-contigs-f.samse-refid";
	string contigs_r_bam_file_name = bam_file_name + "-tea-debug/04-contigs-r.samse-refid";
	string cTea_file_name = bam_file_name + "-tea-debug/05-combined-cTea-refid";

	cout << "Filtering reads without family and combining f and r files" << endl;
	string awk_script_file_name = "/home/el114/kyu/bin/script/teaify-0.3.awk";
	string cmd_awk = awk_script_file_name + " " +
			contigs_f_bam_file_name + " " +
			contigs_r_bam_file_name + " > " + cTea_file_name;
	system(cmd_awk.c_str());

	sort(cTea_file_name);
}

void give_refName(string& bam_file_name) {
	BamReader reader;
	if ( !reader.Open(bam_file_name)) {
		cerr << "Could not open input BAM file" << endl;
		//return false;
	}

	auto& refData = reader.GetReferenceData();
	string line;
	string cTea_refid_file_name = bam_file_name + "-tea-debug/05-combined-cTea-refid-sorted";
	ifstream in_tea(cTea_refid_file_name, ifstream::binary);

	string cTea_file_name = bam_file_name + ".cTea";
	ofstream output(cTea_file_name);

	int64_t first_tab;
	int64_t refID;
	string rest_of_line;

	cout << "Giving refName instead of refID" << endl;
	while(getline(in_tea, line)) {
		first_tab = line.find_first_of("\t");
		if (first_tab == string::npos) {
			continue;
		}
		refID = boost::lexical_cast<int64_t>(line.substr(0, first_tab));
		rest_of_line = line.substr(first_tab);

		output << refData[refID].RefName << rest_of_line << endl;
	}
	output.close();

}


bool ReadClipper::clip_filter_contiggen(
	string& bam_file_name,
	string ref_fa_file_name,
	int64_t minimum_read_length,
	int64_t minimum_base_gap,
	string cap3_options ) {

	mkdir(bam_file_name);

	get_refName(bam_file_name);

	clip(bam_file_name, minimum_read_length);
	run_filter(bam_file_name, minimum_base_gap);

	run_contiggen(bam_file_name, cap3_options, ref_fa_file_name);

	string f_contigs_file_name = bam_file_name + "-tea-debug/04-contigs-f";
	string r_contigs_file_name = bam_file_name + "-tea-debug/04-contigs-r";

	run_bwa(f_contigs_file_name, ref_fa_file_name);
	run_bwa(r_contigs_file_name, ref_fa_file_name);

	filter_family(bam_file_name);

	give_refName(bam_file_name);

	return true;
}

} /* namespace TeaQ */



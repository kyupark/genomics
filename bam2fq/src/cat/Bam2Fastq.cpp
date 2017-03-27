/*
 * Bam2Fastq.cpp
 *
 *  Created on: Mar 24, 2017
 *      Author: el114
 */

#include "Bam2Fastq.hpp"

namespace cat {

Bam2Fastq::Bam2Fastq() :
				n_cores(1) {
	// TODO Auto-generated constructor stub

}

Bam2Fastq::~Bam2Fastq() {
	// TODO Auto-generated destructor stub
}

void Bam2Fastq::process() {
	create_bfi();
	vector<int64_t> block_boundaries;
	load_bfi(block_boundaries);
	write_fastq(block_boundaries);

}

void Bam2Fastq::set_option_parser(const OptionParser& an_option) {
	options = an_option;
	n_cores = options.n_cores;
}

void Bam2Fastq::create_bfi(int64_t n_entries) {

	BamReader reader;
	string bfi_filename = options.input_filename + ".bfi";
	if (boost::filesystem::exists(bfi_filename)) {
		return;
	}
	if (!reader.Open(options.input_filename)) {
		return;
	}

	BamAlignment al;
	int64_t prev_pos=-1;
	int64_t n_local_entries=0;
	auto& bgzf=reader.GetBGZF();
	int64_t n_total_entries=0;
	vector<int64_t> block_boundaries;
	vector<char> dummy;
	while(reader.ReadNextAlignment(dummy)) {
		++n_total_entries;
		++n_local_entries;

		if (n_local_entries > n_entries) {
			string prev_name;
			int64_t cur_pos=-1;

			while(reader.GetNextAlignmentWithName(al)) {

				++n_total_entries;
				++n_local_entries;

				auto& cur_name=al.Name;
				cur_pos=bgzf.Tell();

				if (!prev_name.empty() && prev_name != cur_name) {
					block_boundaries.push_back(prev_pos);
					cout << prev_pos << "\n";
					cout << n_total_entries << " processed\n";
					n_local_entries=0;
					break;
				}
				prev_name=cur_name;
				prev_pos=cur_pos;
			}
		}
	}

	reader.Close();

	block_boundaries.push_back(std::numeric_limits<int64_t>::max());
	ofstream out(bfi_filename, ios::binary);

	for (auto& b : block_boundaries) {
		out << b << "\n";
	}

}
void Bam2Fastq::load_bfi(vector<int64_t>& block_boundaries) {
	string bfi_filename = options.input_filename + ".bfi";
	ifstream in(bfi_filename, ios::binary);

	string line;
	while(getline(in, line)) {
		block_boundaries.push_back(boost::lexical_cast<int64_t>(line));
	}
}

void Bam2Fastq::write_fastq(vector<int64_t>& block_boundaries) {

	// get read groups plus none
	BamReader reader;
	if (!reader.Open(options.input_filename)) {
		return;
	}
	istringstream in(reader.GetHeaderText());
	vector<string> read_groups;
	string line;
	const char* delim = "\t:";
	vector<string> a_cols;
	while (getline(in, line, '\n')) {
		if (string::npos == line.find("@RG")) {
			continue;
		}
		castle::StringUtils::c_string_multi_split(line, delim, a_cols);
		if (a_cols.size() < 2) {
			continue;
		}
		read_groups.push_back(a_cols[2]);
	}
	read_groups.push_back("none");
	reader.Close();

	// Real writing part

//	for (auto& rg : read_groups) {
//		write_read_group(rg, block_boundaries);
//	}

	write_read_group(read_groups[0], block_boundaries);
}

void Bam2Fastq::write_read_group(string& rg, vector<int64_t>& block_boundaries) {
	string output_filename1 = options.output_path + "/" + rg + "_1.fq";
	string output_filename2 = options.output_path + "/" + rg + "_2.fq";

	vector<function<void()> > tasks;
	int64_t n_block_boundaries = block_boundaries.size();
	vector<string> output_files_1(n_block_boundaries - 1);
	vector<string> output_files_2(n_block_boundaries - 1);

	for (int64_t block_id = 0; block_id < n_block_boundaries - 1; ++block_id) {
			tasks.push_back([&, block_id] {
				BamReader local_reader;
				if (!local_reader.Open(options.input_filename)) {
					return;
				}
				auto& local_data = local_reader.GetBGZF();
				int64_t cur_pos = block_boundaries[block_id];
				if(0 != cur_pos) {
					local_data.Seek(cur_pos);
				}
				int64_t next_block_pos = block_boundaries[block_id + 1];

				string str_block_id = boost::lexical_cast<string>(block_id);

				vector<BamAlignment> pair1_alns;
				vector<BamAlignment> pair2_alns;
//				vector<string> pair1_seqs;
//				vector<string> pair2_seqs;
//				vector<string> pair1_quals;
//				vector<string> pair2_quals;

				string prev_read_name;

				BamAlignment al;

				string out_path_1(output_filename1 + "." + str_block_id);
				string out_path_2(output_filename2 + "." + str_block_id);

				ofstream out_disc_1(out_path_1, ios::binary);
				ofstream out_disc_2(out_path_2, ios::binary);
				output_files_1[block_id] = out_path_1;
				output_files_2[block_id] = out_path_2;
//				const bool debug = (0 == block_id);
				const bool debug = false;
				while (local_reader.LoadNextAlignmentCore(al)) {
					auto& cur_name = al.Name;

					if (debug) {
						cout << "[Bam2Fastq.write_read_group] debug-0 " << cur_name << "\n";
					}
					if(prev_read_name != cur_name) {
						if (debug) {
							cout << "[Bam2Fastq.write_read_group] debug-1 " << prev_read_name << "\n";
						}
						if(pair1_alns.size() > 0 && pair2_alns.size() > 0) {
							if (debug) {
								cout << "[Bam2Fastq.write_read_group] debug-2 " << cur_name << "\n";
							}
							// write to files
//							for(uint64_t s_id1 = 0; s_id1 < pair1_seqs.size(); ++s_id1) {
								auto& a_seq = pair1_alns[0].QueryBases;
								auto& a_qual = pair1_alns[0].Qualities;

							string pair1_entry = (boost::format("@%s/1\n%s\n+\n%s\n") % prev_read_name % a_seq % a_qual).str();
							auto& a_seq2 = pair2_alns[0].QueryBases;
							auto& a_qual2 = pair2_alns[0].Qualities;
							string pair2_entry = (boost::format("@%s/2\n%s\n+\n%s\n") % prev_read_name % a_seq2 % a_qual2).str();
							if(debug) {
								cout << "[Bam2Fastq.write_read_group] debug-6 " << a_seq << "\n";
								cout << "[Bam2Fastq.write_read_group] debug-7 " << a_seq2 << "\n";
							}
							out_disc_1 << pair1_entry;
							out_disc_2 << pair2_entry;

//								}
//							}
						}
						prev_read_name = cur_name;
						pair1_alns.clear();
						pair2_alns.clear();
//						pair1_seqs.clear();
//						pair2_seqs.clear();
//						pair1_quals.clear();
//						pair2_quals.clear();
					}
					if("(null)" != cur_name) {
						if (debug) {
							cout << "[Bam2Fastq.write_read_group] debug-3 " << cur_name << "\n";
						}
						if(al.IsFirstMate()) {
							if (debug) {
								cout << "[Bam2Fastq.write_read_group] debug-4 pair 1\n";
							}
							pair1_alns.push_back(al);

						} else if(al.IsSecondMate()) {
							if (debug) {
								cout << "[Bam2Fastq.write_read_group] debug-5 pair 2\n";
							}
							pair2_alns.push_back(al);
						}
					}
					cur_pos = local_data.Tell();
					if(cur_pos >= next_block_pos) {
						break;
					}
				}
				if(pair1_alns.size() > 0 && pair2_alns.size() > 0) {
					auto original_aln1 = pair1_alns;
					auto original_aln2 = pair2_alns;
					// write to files
					auto& a_seq = pair1_alns[0].QueryBases;
					auto& a_qual = pair1_alns[0].Qualities;

					string pair1_entry = (boost::format("@%s/1\n%s\n+\n%s\n") % prev_read_name % a_seq % a_qual).str();
					auto& a_seq2 = pair2_alns[0].QueryBases;
					auto& a_qual2 = pair2_alns[0].Qualities;
					string pair2_entry = (boost::format("@%s/2\n%s\n+\n%s\n") % prev_read_name % a_seq2 % a_qual2).str();
					out_disc_1 << pair1_entry;
					out_disc_2 << pair2_entry;
				}
			});
		}
		castle::ParallelRunner::run_unbalanced_load(n_cores, tasks);
		castle::IOUtils::plain_file_merge(output_filename1, output_files_1, n_cores, false);
		castle::IOUtils::plain_file_merge(output_filename2, output_files_2, n_cores, false);






}


} /* namespace cat */

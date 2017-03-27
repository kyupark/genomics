/*
 * Bam2Fastq.hpp
 *
 *  Created on: Mar 24, 2017
 *      Author: el114
 */

#ifndef CAT_BAM2FASTQ_HPP_
#define CAT_BAM2FASTQ_HPP_

#include "../castle/OptionParser.hpp"
#include "../castle/IOUtils.hpp"
#include "../castle/ParallelRunner.hpp"
#include <api/BamReader.h>
#include <api/BamWriter.h>


namespace cat {

using namespace castle;
using namespace BamTools;

class Bam2Fastq {
public:
	Bam2Fastq();
	~Bam2Fastq();
	void process();
	void set_option_parser(const OptionParser& an_option);
	void create_bfi(int64_t n_entries = 8388608);
	void load_bfi(vector<int64_t>& block_boundaries);
	void write_fastq(vector<int64_t>& block_boundaries);
	void write_read_group(string& rg, vector<int64_t>& block_boundaries);

private:
	OptionParser options;
	int32_t n_cores;

};

} /* namespace cat */

#endif /* CAT_BAM2FASTQ_HPP_ */

/*
 * filter-invisible-positives.hpp
 *
 *  Created on: Feb 24, 2017
 *      Author: el114
 */

#ifndef FILTER_INVISIBLE_POSITIVES_HPP_
#define FILTER_INVISIBLE_POSITIVES_HPP_

#include <iostream>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <api/BamReader.h>


namespace FIP {

using namespace std;
using namespace BamTools;



bool filter_and_save(string& wgs_file_name, string& wes_file_name,
		int64_t range, int64_t minimum_count) {
	BamReader reader;
	if ( !reader.Open(wes_file_name)) {
		cerr << "Could not open input BAM file" << endl;
		//return false;
	}

	ifstream in(wgs_file_name, ifstream::binary);
	string line;
	string chr;
	string pbp;
	string nbp;
	int64_t midpoint=0;
	int64_t alignmentCount=0;
	int64_t refID;
	int64_t leftPosition;
	int64_t rightPosition;

	getline(in, line);
	BamAlignment al;
	while(getline(in, line)){
		if (!line.empty()) {
			stringstream ss(line);
			ss >> chr >> chr >> pbp >> pbp >> pbp >> pbp >> pbp >> nbp;

//			cout << alignmentCount << "\t" << pbp << "\t"  << nbp << "\t" << chr << endl;

			refID = reader.GetReferenceID(chr.substr(3));
			midpoint = ( boost::lexical_cast<int64_t>(pbp) + boost::lexical_cast<int64_t>(nbp) ) / 2;
			leftPosition = midpoint-range;
			rightPosition = midpoint+range;

			reader.OpenIndex(wes_file_name + ".bai");

			if ( reader.HasIndex() ) {
//				cout << "Has Index" << endl;
				// attempt to set region on reader
				if ( !reader.SetRegion(refID, leftPosition, refID, rightPosition) ) {
					cerr << "bamtools count ERROR: set region failed. Check that REGION describes a valid range" << endl;
					reader.Close();
					return false;
				}

				// everything checks out, just iterate through specified region, counting alignments
				while ( reader.GetNextAlignmentCore(al) )
					++alignmentCount;
			}

			// no index data available, we have to iterate through until we
			// find overlapping alignments
			else {
//				cout << "No Index" << endl;
				while ( reader.GetNextAlignmentCore(al) ) {
					if ( (al.RefID == refID)  && ( (al.Position + al.Length) >= leftPosition ) && ( al.Position <= rightPosition) )
					{
						++alignmentCount;
					}
				}
			}

			if (alignmentCount >= minimum_count) {
				cout << alignmentCount << "\t" << line << endl;
			}

			alignmentCount = 0;
		}
	}







	return 0;
}


} /* namespace FIP */

#endif /* FILTER_INVISIBLE_POSITIVES_HPP_ */

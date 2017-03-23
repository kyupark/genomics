// ***************************************************************************
// BamWriter.h (c) 2009 Michael Str�mberg, Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 19 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for producing BAM files
// ***************************************************************************

#ifndef BAMWRITER_H
#define BAMWRITER_H

#include <string>
#include <api/api_global.h>
#include <api/BamAlignment.h>

namespace BamTools {

namespace Internal {
    class BamWriterPrivate;
} // namespace Internal

class API_EXPORT BamWriter {

    // constructor/destructor
    public:
        BamWriter(void);
        ~BamWriter(void);

    // public interface
    public:
        // closes the alignment archive
        void Close(void);
        // opens the alignment archive
        bool Open(const std::string& filename, 
                  const std::string& samHeader, 
                  const BamTools::RefVector& referenceSequences, 
                  bool writeUncompressed = false);
        // opens the SAM format file
		bool SAMOpen(const std::string& filename,
				  const std::string& samHeader,
				  const BamTools::RefVector& referenceSequences);
		bool SAMOpenNoHeader(const std::string& filename,
						  const BamTools::RefVector& referenceSequences);
        // saves the alignment to the alignment archive
        void SaveAlignment(const BamTools::BamAlignment& al);
        void SaveAlignmentBigEndianCoreOnly(const BamTools::BamAlignment& al);
		void SaveAlignmentLittleEndianCoreOnly(const BamTools::BamAlignment& al);
		void SaveAlignmentBigEndianNonCore(const BamTools::BamAlignment& al);
		void SaveAlignmentLittleEndianNonCore(const BamTools::BamAlignment& al);
        void SaveSAMAlignment(const BamTools::BamAlignment& al);
        static std::string GetSAMAlignment(const BamTools::BamAlignment& a, const RefVector& m_references);

    // private implementation
    private:
	Internal::BamWriterPrivate* d;
};

} // namespace BamTools

#endif // BAMWRITER_H

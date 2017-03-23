// ***************************************************************************
// BamWriter.cpp (c) 2009 Michael Str�mberg, Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 22 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for producing BAM files
// ***************************************************************************

#include "../../src/api/BamWriter.h"

#include "../../src/api/internal/BamWriter_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <iostream>
using namespace std;

// constructor
BamWriter::BamWriter(void) {
    d = new BamWriterPrivate;
}

// destructor
BamWriter::~BamWriter(void) {
    delete d;
    d = 0;
}

// closes the alignment archive
void BamWriter::Close(void) {
  d->Close();
}

// opens the alignment archive
bool BamWriter::Open(const string& filename,
		     const string& samHeader,
		     const RefVector& referenceSequences,
		     bool isWriteUncompressed)
{
    return d->Open(filename, samHeader, referenceSequences, isWriteUncompressed);
}

bool BamWriter::SAMOpen(const std::string& filename,
				  const std::string& samHeader,
				  const BamTools::RefVector& referenceSequences) {
	return d->SAMOpen(filename, samHeader, referenceSequences);
}
bool BamWriter::SAMOpenNoHeader(const std::string& filename,
						  const BamTools::RefVector& referenceSequences) {
	return d->SAMOpenNoHeader(filename, referenceSequences);
}
// saves the alignment to the alignment archive
void BamWriter::SaveAlignment(const BamAlignment& al) {
    d->SaveAlignment(al);
}
void BamWriter::SaveAlignmentBigEndianCoreOnly(const BamTools::BamAlignment& al) {
	d->SaveAlignmentBigEndianCoreOnly(al);
}
void BamWriter::SaveAlignmentLittleEndianCoreOnly(const BamTools::BamAlignment& al) {
	d->SaveAlignmentLittleEndianCoreOnly(al);
}
void BamWriter::SaveAlignmentBigEndianNonCore(const BamTools::BamAlignment& al) {
	d->SaveAlignmentBigEndianNonCore(al);
}
void BamWriter::SaveAlignmentLittleEndianNonCore(const BamTools::BamAlignment& al) {
	d->SaveAlignmentLittleEndianNonCore(al);
}
void BamWriter::SaveSAMAlignment(const BamTools::BamAlignment& al) {
	d->SaveSAMAlignment(al);
}

string BamWriter::GetSAMAlignment(const BamTools::BamAlignment& al, const RefVector& m_references) {
	return Internal::BamWriterPrivate::GetSAMAlignment(al, m_references);
}

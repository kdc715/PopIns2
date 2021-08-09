#include <sstream>
#include <cerrno>

#include <seqan/file.h>
#include <seqan/sequence.h>

#include "util.h"
#include "argument_parsing.h"
#include "crop_unmapped.h"

#ifndef POPINS2_MERGE_SETMATE_H_
#define POPINS2_MERGE_SETMATE_H_

using namespace seqan;



// ==========================================================================

inline void setMates(BamAlignmentRecord & record1, BamAlignmentRecord & record2) {
    SEQAN_ASSERT(!hasFlagFirst(record1) || !hasFlagFirst(record2));
    SEQAN_ASSERT(!hasFlagLast(record1) || !hasFlagLast(record2));

    // Set the next ref id.
    record1.rNextId = record2.rID;
    record2.rNextId = record1.rID;

    // Set the next ref pos.
    record1.pNext = record2.beginPos;
    record2.pNext = record1.beginPos;

    // Fix the next unmapped flag.
    if (hasFlagUnmapped(record2)) record1.flag |= BAM_FLAG_NEXT_UNMAPPED;
    else record1.flag &= ~BAM_FLAG_NEXT_UNMAPPED;
    if (hasFlagUnmapped(record1)) record2.flag |= BAM_FLAG_NEXT_UNMAPPED;
    else record2.flag &= ~BAM_FLAG_NEXT_UNMAPPED;

    // Fix the next reversed flag.
    if (hasFlagRC(record2)) record1.flag |= BAM_FLAG_NEXT_RC;
    else record1.flag &= ~BAM_FLAG_NEXT_RC;
    if (hasFlagRC(record1)) record2.flag |= BAM_FLAG_NEXT_RC;
    else record2.flag &= ~BAM_FLAG_NEXT_RC;

    // Fix first/second in pair flags.
    if (hasFlagFirst(record1)) record2.flag |= BAM_FLAG_LAST;
    if (hasFlagFirst(record2)) record1.flag |= BAM_FLAG_LAST;
    if (hasFlagLast(record1)) record2.flag |= BAM_FLAG_FIRST;
    if (hasFlagLast(record2)) record1.flag |= BAM_FLAG_FIRST;

    // Set flag paired.
    record1.flag |= BAM_FLAG_MULTIPLE;
    record2.flag |= BAM_FLAG_MULTIPLE;
}

// ==========================================================================

// Correct the reference ids of a BamAlignmentRecord for the concatenated header.
template<typename TNameStore>
inline void readRecordAndCorrectRIds(BamAlignmentRecord & record, BamFileIn & stream, NameStoreCache<TNameStore> & nameStoreCache){
    readRecord(record, stream);

    if (record.rID != BamAlignmentRecord::INVALID_REFID)
    {
        CharString rName = contigNames(context(stream))[record.rID];
        getIdByName(record.rID, nameStoreCache, rName);
    }
    if (record.rNextId != BamAlignmentRecord::INVALID_REFID)
    {
        CharString rNextName = contigNames(context(stream))[record.rNextId];
        getIdByName(record.rNextId, nameStoreCache, rNextName);
    }
}

// ==========================================================================

inline void mergeHeaders(BamHeader & header, FormattedFileContext<BamFileOut, Owner<> >::Type & context, BamFileIn & stream1, BamFileIn & stream2){
    StringSet<CharString> contigNames;
    NameStoreCache<StringSet<CharString> > nameStoreCache;
    String<int32_t> contigLengths;

    // Read and append the two headers. Remove duplicate entries.
    readHeader(header, stream1);
    BamHeader header2;
    readHeader(header2, stream2);
    for (unsigned i = 0; i < length(header2); ++i)
    {
        if (header2[i].type != BAM_HEADER_FIRST)
            appendValue(header, header2[i]);
    }
    std::stable_sort(begin(header, Standard()), end(header, Standard()), BamHeaderRecordTypeLess());

    // Fill sequence names into nameStoreCache.
    for (unsigned i = 0; i < length(header); ++i)
    {
        if (header[i].type == BAM_HEADER_REFERENCE)
        {
            CharString name, len;
            for (unsigned j = 0; j < length(header[i].tags); ++j)
            {
                if (header[i].tags[j].i1 == "SN")
                    name = header[i].tags[j].i2;
                else if (header[i].tags[j].i1 == "LN")
                    len = header[i].tags[j].i2;
            }
            appendName(context._contigNamesCache, name);
            int32_t l;
            lexicalCast<int32_t>(l, len);
            appendValue(context._contigLengths, l);
        }
    }
}

// ==========================================================================

// This function is adapted from samtools code (fuction strnum_cmp in bam_sort.c) to ensure the exact same sort order.
int
compare_qName(CharString & nameA, CharString & nameB)
{
    const char * _a = toCString(nameA);
    const char * _b = toCString(nameB);
    const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
    const unsigned char *pa = a, *pb = b;
    while (*pa && *pb) {
        if (isdigit(*pa) && isdigit(*pb)) {
            while (*pa == '0') ++pa;
            while (*pb == '0') ++pb;
            while (isdigit(*pa) && isdigit(*pb) && *pa == *pb) ++pa, ++pb;
            if (isdigit(*pa) && isdigit(*pb)) {
                int i = 0;
                while (isdigit(pa[i]) && isdigit(pb[i])) ++i;
                return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
            } else if (isdigit(*pa)) return 1;
            else if (isdigit(*pb)) return -1;
            else if (pa - a != pb - b) return pa - a < pb - b? 1 : -1;
        } else {
            if (*pa != *pb) return (int)*pa - (int)*pb;
            ++pa; ++pb;
        }
    }
    return *pa? 1 : *pb? -1 : 0;
}

// ==========================================================================
// Function merge_and_set_mate()
// ==========================================================================

bool
merge_and_set_mate(CharString &mergedBam, unsigned &nonContigSeqs, CharString &nonRefBam, CharString &remappedBam)
{
    std::ostringstream msg;
    msg << "Merging bam files " << nonRefBam << " and " << remappedBam;
    printStatus(msg);

    // Open the two input streams (can read SAM and BAM files).
    BamFileIn nonRefStream(toCString(nonRefBam));
    BamFileIn remappedStream(toCString(remappedBam));

    printStatus(" - merging headers...");

    // Prepare a header for the output file.
    BamHeader outHeader;
    FormattedFileContext<BamFileOut, Owner<> >::Type bamContext;
    mergeHeaders(outHeader, bamContext, nonRefStream, remappedStream);

    printStatus(" - writing header...");

    // Open the output stream and write the header.
    FormattedFileContext<BamFileOut, Dependent<> >::Type bamContextDep(bamContext);
    BamFileOut outStream(bamContextDep, toCString(mergedBam));
    writeHeader(outStream, outHeader);

    nonContigSeqs = length(contigNames(context(nonRefStream)));

    printStatus(" - merging read records...");

    // Read the first record from each input file. Correct ids in records from remappedStreams for new header.
    BamAlignmentRecord record1, record2;
    if (!atEnd(nonRefStream)) readRecordAndCorrectRIds(record1, nonRefStream, contigNamesCache(bamContextDep));
    else record1.qName = "*";
    if (!atEnd(remappedStream)) readRecordAndCorrectRIds(record2, remappedStream, contigNamesCache(bamContextDep));
    else record2.qName = "*";

    // Iterate both input files, set mate positions in pairs, and write all records to the output file.
    while (record1.qName != "*" || record2.qName != "*")
    {
        while ((compare_qName(record2.qName, record1.qName) < 0 || record1.qName == "*") && record2.qName != "*")
        {
            writeRecord(outStream, record2);
            if (!atEnd(remappedStream)) readRecordAndCorrectRIds(record2, remappedStream, contigNamesCache(bamContextDep));
            else record2.qName = "*";
        }

        bool incr1 = false;
        while (record1.qName == record2.qName && record2.qName != "*")
        {
            incr1 = true;
            setMates(record1, record2);
            writeRecord(outStream, record1);
            writeRecord(outStream, record2);
            if (!atEnd(remappedStream)) readRecordAndCorrectRIds(record2, remappedStream, contigNamesCache(bamContextDep));
            else record2.qName = "*";
        }
        if (incr1)
        {
            if (!atEnd(nonRefStream)) readRecordAndCorrectRIds(record1, nonRefStream, contigNamesCache(bamContextDep));
            else record1.qName = "*";
        }

        while ((compare_qName(record1.qName, record2.qName) < 0 || record2.qName == "*") && record1.qName != "*")
        {
            writeRecord(outStream, record1);
            if (!atEnd(nonRefStream)) readRecordAndCorrectRIds(record1, nonRefStream, contigNamesCache(bamContextDep));
            else record1.qName = "*";
        }
    }

    return 0;
}

// ==========================================================================
// Function popins2_merge_and_set_mate()
// ==========================================================================

bool
popins2_merge_and_set_mate(int argc, char const ** argv)
{
    // Parse the command line to get option values.
    MergeSetMateOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    CharString workingDirectory = getFileName(options.prefix, options.sampleID);
    
    CharString mergedBam = getFileName(workingDirectory, "merged.bam");
    unsigned nonContigSeqs = 0;
    CharString nonRefBam = getFileName(workingDirectory, "non_ref.bam");
    CharString remappedBam = getFileName(workingDirectory, "remapped.bam");

    bool ret = merge_and_set_mate(mergedBam, nonContigSeqs, nonRefBam, remappedBam);
    return ret;
}


#endif // #ifndef POPINS2_MERGE_SETMATE_H_

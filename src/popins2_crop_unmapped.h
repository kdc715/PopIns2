#include <sstream>
#include <cerrno>

#include <seqan/file.h>
#include <seqan/sequence.h>

#include "util.h"
#include "argument_parsing.h"
#include "crop_unmapped.h"

#ifndef POPINS2_CROP_UNMAPPED_H_
#define POPINS2_CROP_UNMAPPED_H_

using namespace seqan;

// ==========================================================================


bool retrieveSampleID(CharString & sampleID, CharString & mappingBam)
{
    BamFileIn inStream(toCString(mappingBam));

    BamHeader header;
    readHeader(header, inStream);

    for (unsigned i = 0; i < length(header); ++i)
    {
        if (header[i].type != BamHeaderRecordType::BAM_HEADER_READ_GROUP)
            continue;

        for (unsigned j = 0; j < length(header[i].tags); ++j)
        {
            if (header[i].tags[j].i1 != "SM")
                continue;

            sampleID = header[i].tags[j].i2;
            return 0;
        }
    }
    std::cerr << "ERROR: Could not find sample ID in BAM file header." << std::endl;
    return 1;
}


// ==========================================================================
// Function popins2_crop_unmapped()
// ==========================================================================

int popins2_crop_unmapped(int argc, char const ** argv)
{
    std::ostringstream msg;

    // Parse the command line to get option values.
    CropUnmappedOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Retrieve the sample ID from the first read group listed in BAM file header.
    if (options.sampleID == "" && retrieveSampleID(options.sampleID, options.mappingFile) == 1)
        return 7;

    // Create working directory if it does not exist.
    CharString workingDirectory = getFileName(options.prefix, options.sampleID);
    if (mkdir(toCString(workingDirectory), 0755) == 0)
    {
        msg.str("");
        msg << "Working directory created at " << workingDirectory;
        printStatus(msg);
    }

    SampleInfo info = initSampleInfo(options.mappingFile, options.sampleID, options.adapters);

    float as_factor = options.alignment_score_factor;

    CharString matesBam = getFileName(workingDirectory, "mates.bam");

    CharString fastqFirst = getFileName(workingDirectory, "paired.1.fastq");
    CharString fastqSecond = getFileName(workingDirectory, "paired.2.fastq");
    CharString fastqSingle = getFileName(workingDirectory, "single.fastq");
    Triple<CharString> fastqFiles = Triple<CharString>(fastqFirst, fastqSecond, fastqSingle);



    // check if files already exits
    std::fstream stream(toCString(fastqFirst));
    if (!stream.is_open())
    {
        msg.str("");
        msg << "Cropping unmapped reads from " << options.mappingFile;
        printStatus(msg);

        // Crop unmapped reads and reads with unreliable mappings from the input bam file.
        if (options.adapters == "HiSeqX")
        {
            if (crop_unmapped(info.avg_cov, fastqFiles, matesBam, options.mappingFile, options.humanSeqs, HiSeqXAdapters(), as_factor) != 0)
                return 7;
        }
        else if (options.adapters == "HiSeq")
        {
            if (crop_unmapped(info.avg_cov, fastqFiles, matesBam, options.mappingFile, options.humanSeqs, HiSeqAdapters(), as_factor) != 0)
                return 7;
        }
        else
        {
            if (crop_unmapped(info.avg_cov, fastqFiles, matesBam, options.mappingFile, options.humanSeqs, NoAdapters(), as_factor) != 0)
                return 7;
        }

        CharString sampleInfoFile = getFileName(workingDirectory, "POPINS_SAMPLE_INFO");
        writeSampleInfo(info, sampleInfoFile);

        msg.str("");
        msg << "Sample info written to \'" << sampleInfoFile << "\'.";
        printStatus(msg);

      
    }

    return res;

}

#endif // #ifndef POPINS2_CROP_UNMAPPED_H_


#include <sstream>
#include <cerrno>

#include <seqan/file.h>
#include <seqan/sequence.h>

#include "util.h"
#include "argument_parsing.h"
#include "crop_unmapped.h"

#ifndef POPINS2_REMAPPING_H_
#define POPINS2_REMAPPING_H_

using namespace seqan;


// ==========================================================================
// Function popins2_remapping()
// could be taken apart in snakemake, only calling external functions(SAMTOOLS/BWA)
// ==========================================================================

inline int popins2_remapping(int argc, char const ** argv)
{
    std::stringstream cmd;
    std::ostringstream msg;

    RemappingOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
    {
        return res;
    } 
    // Create working directory if it does not exist.
    CharString workingDirectory = getFileName(options.prefix, options.sampleID);
    if (mkdir(toCString(workingDirectory), 0755) == 0)
    {
        msg.str("");
        msg << "Working directory created at " << workingDirectory;
        printStatus(msg);
    }

    CharString fastqFirst = getFileName(workingDirectory, "paired.1.fastq");
    CharString fastqSecond = getFileName(workingDirectory, "paired.2.fastq");
    CharString fastqSingle = getFileName(workingDirectory, "single.fastq");
    Triple<CharString> fastqFiles = Triple<CharString>(fastqFirst, fastqSecond, fastqSingle);
    
    Triple<CharString> fastqFilesTemp = fastqFiles;
       

    CharString f1 = options.prefix;
    f1 += "remapped.sam";
    CharString remappedSam = getFileName(options.workingDir, f1);

    CharString f2 = options.prefix;
    f2 += "remapped.bam";
    CharString remappedBam = getFileName(options.workingDir, f2);

    CharString f3 = options.prefix;
    f3 += "remapped.bam.bai";
    CharString remappedBai = getFileName(options.workingDir, f3);

    CharString f4 = options.prefix;
    f4 += "remapped_unsorted.bam";
    CharString remappedUnsortedBam = getFileName(options.workingDir, f4);

    
    msg << "Remapping unmapped reads using " << BWA;
    printStatus(msg);

    

    // Run BWA on unmapped reads (pairs).
    cmd.str("");
    cmd << BWA << " mem -t " << options.threads << " " << options.referenceFile << " " << fastqFilesTemp.i1 << " " << fastqFilesTemp.i2 << " > " << remappedSam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while running bwa on " << fastqFilesTemp.i1 << " and " << fastqFilesTemp.i2 << std::endl;
        return 1;
    }

    remove(toCString(fastqFilesTemp.i1));
    remove(toCString(fastqFilesTemp.i2));

    // Run BWA on unmapped reads (single end).
    cmd.str("");
    cmd << BWA << " mem -t " << options.threads << " " << options.referenceFile << " " << fastqFilesTemp.i3 << " | awk '$1 !~ /@/' >> " << remappedSam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while running bwa on " << fastqFilesTemp.i3 << std::endl;
        return 1;
    }

    remove(toCString(fastqFilesTemp.i3));

    msg.str("");
    msg << "Converting BWA output " << remappedSam << " to bam format.";
    printStatus(msg);

    // Convert BWA output to bam.
    cmd.str("");
    cmd << SAMTOOLS << " view -S -h -b " << remappedSam << " > " << remappedUnsortedBam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while converting BWA output " << remappedSam << " to bam." << std::endl;
        return 1;
    }
    remove(toCString(remappedSam));

    msg.str("");
    msg << "Sorting " << remappedUnsortedBam << " using " << SAMTOOLS;
    printStatus(msg);

    // Sort bam file.
    cmd.str("");
    cmd << SAMTOOLS << " sort -@ " << options.threads << " -m " << options.memory << " -o " << remappedBam << " " << remappedUnsortedBam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while sorting BWA output " << remappedUnsortedBam << std::endl;
        return 1;
    }

    msg.str("");
    msg << "Indexing " << remappedBam << " using " << SAMTOOLS;
    printStatus(msg);

    // Index bam file.
    cmd.str("");
    cmd << SAMTOOLS << " index " << remappedBam;
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while indexing BWA output " << remappedBam << std::endl;
        return 1;
    }

    msg.str("");
    msg << "Cropping unmapped reads from " << remappedBam;
    printStatus(msg);

    // Crop unmapped and create bam file of remapping.
    if (crop_unmapped(fastqFiles, remappedUnsortedBam, remappedBam, options.humanSeqs, NoAdapters(), options.alignment_score_factor) != 0)
        return 1;
    remove(toCString(remappedBai));

    msg.str("");
    msg << "Sorting " << remappedUnsortedBam << " by read name using " << SAMTOOLS;
    printStatus(msg);

    // Sort <WD>/remapped.bam by read name.
    cmd.str("");
    cmd << SAMTOOLS << " sort -n -@ " << options.threads << " -m " << options.memory << " -o " << remappedBam << " " << remappedUnsortedBam << " ";
    if (system(cmd.str().c_str()) != 0)
    {
        std::cerr << "ERROR while sorting " << remappedUnsortedBam << std::endl;
        return 1;
    }

    remove(toCString(remappedUnsortedBam));

    return 0;
}



#endif // #ifndef POPINS2_REMAPPING_H_


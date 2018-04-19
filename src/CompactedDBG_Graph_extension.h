/*!
* \file    src/CompactedDBG_Graph_extension.h
* \brief   Library for a compacted de Bruijn Graph using data extension.
*
*/
#ifndef COMPACTED_DBG_GRAPH_EXTENSION_
#define COMPACTED_DBG_GRAPH_EXTENSION_

#include "CompactedDBG_Data_extension.h"

#include <cstdint>                      /* uint8_t */
#include <seqan/seq_io.h>

#include "../../prettyprint/prettyprint.h"      // delete for release



// =========================
// Global
// =========================
typedef std::vector<std::vector<unsigned int> > BubblePathSet;
typedef std::vector<unsigned int> BubblePath;


// =========================
// Struct
// =========================
/*!
* \class        ExtendedCDBG
* \headerfile   src/CompactedDBG_Graph_extension.h
* \brief        Struct to store a compacted DBG plus data extension.
*/
struct ExtendedCDBG : public CompactedDBG<DataExtension> {

    public:
        // constructor
        ExtendedCDBG(int kmer_length, int minimizer_length) :   CompactedDBG<DataExtension,void> (kmer_length, minimizer_length),
                                                                DFS_MAX_DIST(kmer_length<<1),
                                                                id_init_status(false)
                                                                {}

        void init_ids();

        bool annotate_kmer_coverage(const vector<string> &sample_fastx_names);

        bool small_bubble_removal();

    private:

        uint8_t whereToGo(const UnitigMap<DataExtension> &um, const UnitigMap<DataExtension> &src) const;

        uint8_t getDegree(const UnitigMap<DataExtension> &um) const;  // return overall degree (in- plus outdegree) of a unitig

        void clear_traversal_attributes();

        bool DFS_Direction_Init(const UnitigMap<DataExtension> &um, const uint8_t direction);
        bool DFS_Direction_Recursion(const UnitigMap<DataExtension> &um, const UnitigMap<DataExtension> &src, unsigned int dist);

        bool Traceback_Init(const UnitigMap<DataExtension> &src, const UnitigMap<DataExtension> &traceback_src);
        bool Traceback_Visit(const UnitigMap<DataExtension> &um, const UnitigMap<DataExtension> &src, BubblePath &path, uint8_t &nb_branchingBubblePaths);

        bool delete_bubble_path(const BubblePathSet &bps, const uint8_t nb_branchingBubblePaths);

        const static uint8_t GO_FORWARD = 0x0;
        const static uint8_t GO_BACKWARD = 0x1;

        const size_t DFS_MAX_DIST;   // default: two times k-mer length

        bool id_init_status;

};





















#endif /*COMPACTED_DBG_GRAPH_EXTENSION_*/
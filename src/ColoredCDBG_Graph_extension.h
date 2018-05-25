/*!
* \file    src/ColoredCDBG_Graph_extension.h
* \brief   Library for a colored compacted de Bruijn Graph using unitig extension.
*
*/
#ifndef COLORED_CDBG_GRAPH_EXTENSION_
#define COLORED_CDBG_GRAPH_EXTENSION_


#include <seqan/seq_io.h>
#include <seqan/misc/union_find.h>

#include "ColoredCDBG_Data_extension.h"

typedef std::vector<std::vector<UnitigColorMap<UnitigExtension> > > PathSet;
typedef std::vector<UnitigColorMap<UnitigExtension> > UnitigPath;



// =========================
// Struct
// =========================
/*!
* \class        ExtendedCCDBG
* \headerfile   src/ColoredCDBG_Graph_extension.h
* \brief        Struct to store a colored compacted DBG plus unitig extensions.
*/
struct ExtendedCCDBG : public ColoredCDBG<UnitigExtension> {

    private:

        bool id_init_status;

        seqan::UnionFind<unsigned> UF;

        const static uint8_t GO_FORWARD = 0x0;
        const static uint8_t GO_BACKWARD = 0x1;

        bool DFS_CLEAN = true;

    public:
        ExtendedCCDBG(int kmer_length = 31, int minimizer_length = 23);      // hidden inits! (see definition)

        void init_ids();
        void print_ids();
        bool is_id_init() const {return id_init_status;}

        bool connected_components(const CCDBG_Build_opt &graph_options);
        size_t count_connected_components();
        seqan::UnionFind<unsigned> getUF() const {return UF;}

        float entropy(const std::string &sequence);

        bool is_DFS_clean() const {return DFS_CLEAN;}
        void make_DFS_clean();

        uint8_t whereToGo(const UnitigColorMap<UnitigExtension> &um, const UnitigColorMap<UnitigExtension> &src) const;

        bool DFS_Init(const UnitigColorMap<UnitigExtension> &ucm, const uint8_t direction, const bool verbose);
        bool DFS_Visit(const UnitigColorMap<UnitigExtension> &ucm,
                       const UnitigColorMap<UnitigExtension> &src,
                       const UnitigColorMap<UnitigExtension> &anchor,
                       const bool verbose);

};












#endif /*COLORED_CDBG_GRAPH_EXTENSION_*/
#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <bifrost/CompactedDBG.hpp>

#include <unordered_map>

#include "../src/ExtendedCDBG.h"
#include "../src/argument_parsing.h"



CDBG_Build_opt graph_opt;
SEQAN_DEFINE_TEST(test_bifrost_parameter){
    // -----------------------------
    // | Test graph build options  |
    // -----------------------------
    std::vector<std::string> infiles;
    infiles.push_back("./testdata/S0001_human_simulated.fq");
    infiles.push_back("./testdata/S0002_human_simulated.fq");
    graph_opt.filename_in = infiles;
    graph_opt.prefixFilenameOut = "union_test_out";
    graph_opt.nb_unique_kmers = 2488350;
    graph_opt.nb_non_unique_kmers = 874940;
    graph_opt.nb_threads = 4;
    graph_opt.clipTips = true;
    graph_opt.deleteIsolated = true;
    SEQAN_ASSERT_EQ(check_ProgramOptions(graph_opt), true);
}

ExtendedCDBG g(graph_opt.k, graph_opt.g);
SEQAN_DEFINE_TEST(test_bifrost_graphfunctions){
    // -----------------------
    // | Test Bifrost methods  |
    // -----------------------
    g.build(graph_opt);
    g.simplify(graph_opt.deleteIsolated, graph_opt.clipTips, graph_opt.verbose);
    SEQAN_ASSERT_EQ(g.size(), 119u);
    SEQAN_ASSERT_EQ(g.write(graph_opt.prefixFilenameOut, graph_opt.nb_threads, true, graph_opt.verbose), true);
}

SEQAN_DEFINE_TEST(test_init_ids){
    g.init_ids();
    SEQAN_ASSERT_EQ(g.is_init(), true);
}

SEQAN_DEFINE_TEST(test_connectedcomponents){
    bool cc_build = g.connected_components(graph_opt);
    SEQAN_ASSERT_EQ(cc_build, true);
    size_t nb_cc = g.count_connected_components();
    SEQAN_ASSERT_EQ(nb_cc, 93u);
}

SEQAN_DEFINE_TEST(test_neighbors_and_bit_operations){
    size_t nb_potential_splitnodes = 0;
    for (auto &unitig : g){
        // -----------------------------
        // | Test amount of neighbors  |
        // -----------------------------
        ForwardCDBG<UnitigExtension, false> fw_dbg = unitig.getSuccessors();
        BackwardCDBG<UnitigExtension, false> bw_dbg = unitig.getPredecessors();

        size_t i=0;
        for (auto &suc : fw_dbg) ++i;

        size_t j=0;
        for (auto &pre : bw_dbg) ++j;

        if (i>1 && j>1){
            SEQAN_ASSERT_EQ(i, 2u);
            SEQAN_ASSERT_EQ(j, 2u);
            ++nb_potential_splitnodes;
        }
        // ------------------------
        // | Test base accession  |
        // ------------------------
        if (unitig.getData()->getID() == 4) {
            UnitigMap<UnitigExtension>::neighbor_iterator suc1 = fw_dbg.begin();
            UnitigMap<UnitigExtension>::neighbor_iterator pre1 = bw_dbg.begin();

            char firstPredecessor_lastBase = (pre1->getTail()).getChar(graph_opt.k-1);   // getChar(offset) returns the base from Kmer at offset from the beginning
            char firstSuccessor_firstBase = (suc1->getHead()).getChar(0);                   // getChar(offset) returns the base from Kmer at offset from the beginning

            SEQAN_ASSERT_EQ(firstPredecessor_lastBase, 'T');
            SEQAN_ASSERT_EQ(firstSuccessor_firstBase, 'C');

            // -------------------------------------------
            // | Test succint bit storage and bit masks  |
            // -------------------------------------------
            uint8_t pre_bitmask = bitmask_encoder_predecessor[firstPredecessor_lastBase];
            uint8_t suc_bitmask = bitmask_encoder_successor[firstSuccessor_firstBase];
            uint8_t test_bit_encoding = setNeighborPairFromBases(pre_bitmask, suc_bitmask);
            SEQAN_ASSERT_EQ(test_bit_encoding, 0b00001101);

            pair<uint8_t, uint8_t> test_bit_decoding = getBasesFromNeighborPair(test_bit_encoding);
            char f = bitmask_decoder[test_bit_decoding.first];
            char s = bitmask_decoder[test_bit_decoding.second];
            SEQAN_ASSERT_EQ(f, 'T');
            SEQAN_ASSERT_EQ(s, 'C');
        }
    }
    SEQAN_ASSERT_EQ(nb_potential_splitnodes, 5u);
}

/* NOTE:    This version of storing pointer doesn't work.
 *          This test is currently not called in the the testsuite.
 *          Instead store copies of iterator (see test_pointer_map2).
 */
SEQAN_DEFINE_TEST(test_pointer_map){
    // build phase
    std::unordered_map<unsigned, UnitigMap<UnitigExtension>*> m;
    CompactedDBG<UnitigExtension>::iterator it = g.begin(), it_end = g.end();
    for (; it != it_end; ++it){
        unsigned id = it->getData()->getID();
        UnitigMap<UnitigExtension>* p_um = &(*it);
        cout << "ID: " << p_um->getData()->getID() << " | Pointer: " << p_um << endl;
        m[id] = p_um;
    }
    SEQAN_ASSERT_EQ(m.size(), 119u);

    // access phase
    UnitigMap<UnitigExtension>* test_p_um = m[3];
    unsigned test_id = test_p_um->getData()->getID();
    SEQAN_ASSERT_EQ(test_id, 3u);   // CRASHES!!!
}

/* This is a test whether the concept of test_pointer_map works in general.
 */
SEQAN_DEFINE_TEST(test_pointer_map2){
    // build phase
    seqan::String<unsigned> testString;
    seqan::resize(testString, 10, 0);
    unsigned i = 0;
    for (seqan::Iterator<seqan::String<unsigned> >::Type it = seqan::begin(testString); it != seqan::end(testString); ++it)
    {
        *it = i;
        ++i;
    }

    // access phase
    seqan::Iterator<seqan::String<unsigned> >::Type it = seqan::begin(testString), it_end = seqan::end(testString);
    for (; it != it_end; ++it){
        unsigned * p = &*it;
        cout << "ItValue=" << *it << "\titAdr=" << it << "\tPtrVal=" << *p << "\tPtrAdr=" << p << endl;
    }

}

SEQAN_DEFINE_TEST(test_pointer_map3){
    // build phase
    std::unordered_map<unsigned, CompactedDBG<UnitigExtension>::iterator> m;
    CompactedDBG<UnitigExtension>::iterator it = g.begin(), it_end = g.end();
    for (; it != it_end; ++it){
        unsigned id = it->getData()->getID();
        CompactedDBG<UnitigExtension>::iterator p_um = it;
        //cout << "ID: " << p_um->getData()->getID() << " | Pointer: " << &(*p_um) << endl;
        m[id] = p_um;
    }
    SEQAN_ASSERT_EQ(m.size(), 119u);

    // access phase
    CompactedDBG<UnitigExtension>::iterator test_i_um = m[4];
    unsigned test_id = test_i_um->getData()->getID();
    SEQAN_ASSERT_EQ(test_id, 4u);
    SEQAN_ASSERT_EQ(test_i_um->toString(), "CCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGA");
}



SEQAN_BEGIN_TESTSUITE(test_popins2){
	// call tests here
    SEQAN_CALL_TEST(test_bifrost_parameter);
    SEQAN_CALL_TEST(test_bifrost_graphfunctions);
    SEQAN_CALL_TEST(test_init_ids);
    SEQAN_CALL_TEST(test_connectedcomponents);
    SEQAN_CALL_TEST(test_neighbors_and_bit_operations);
    SEQAN_CALL_TEST(test_pointer_map3);
}
SEQAN_END_TESTSUITE

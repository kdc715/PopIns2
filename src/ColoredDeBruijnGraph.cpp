#include "ColoredDeBruijnGraph.h"

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <queue>
#include <fstream>
#include <map>


uint8_t Traceback::recursion_priority_counter = 0;



inline void Traceback::join(const Traceback &t){

    for (auto it = t.cbegin(); it != t.cend(); ++it)
        this->push_back(*it);
}


inline void Traceback::cutconcat(string &s, const VSequences &path, const size_t k) const{
    // since the traceback stored the sequences from sink to source, the concatination has to run backwards over the vector
    for (auto rit = path.crbegin(); rit != path.crend(); ++rit){
        if (rit == path.crbegin())
            s+=(*rit);
        else
            s+=rit->substr(k-1);
    }
}


/*!
 * \fn      bool Traceback::write(const std::string &filename, ofstream &ofs) const
 * \return  true if write was successful
 */
bool Traceback::write(ofstream &ofs, const size_t k, size_t &counter) const{
    if (ofs.is_open()){
        // loop through all traceback paths found in one node
        for (auto it=cbegin(); it!=cend(); ++it){
            std::string seq;
            cutconcat(seq, *it, k);
            ++counter;

            ofs << ">contig_" << counter << "\n";
            ofs << seq << "\n";
        }
        return 1;
    }
    else
        cout << "Error: Unable to open file.";
    return 0;
}


// default constructor
ExtendedCCDBG::ExtendedCCDBG(int kmer_length, int minimizer_length) :   ColoredCDBG<UnitigExtension> (kmer_length, minimizer_length),
                                                                        id_init_status(false) {
    /* 1) IDs are not initiated at construction time (see init_ids())
     * 2) The UnionFind vector will be empty an construction time, will be resized at use.
     */
}


void ExtendedCCDBG::init_ids(){
    size_t i=1;           // starting index is 1 because that's how Bifrost counts
    for (auto &unitig : *this){
        DataAccessor<UnitigExtension>* da = unitig.getData();
        UnitigExtension* ue = da->getData(unitig);      // ue is a POINTER to a UnitigExtension
        ue->setID(i);
        ++i;
    }
    id_init_status = true;
}


void ExtendedCCDBG::print_ids(){
    if (is_id_init() == true){
        std::cout << "[PRINT] ";
        for (auto &unitig : *this){
            DataAccessor<UnitigExtension>* da = unitig.getData();
            UnitigExtension* ue = da->getData(unitig);
            unsigned id = ue->getID();
            std::cout << id << ", ";
        }
        std::cout << std::endl;
    }
    else
        cerr << "[WARNING] Unitig IDs were not printed because they are not initialized." << endl;
}


/*!
 * \fn      size_t ExtendedCCDBG::count_connected_components()
 * \brief   This function computes the amount of distict connected components in the graph. The connected_components()
 *          function need to be ran successfully to call count_connected_components(). This function is mainly for
 *          debug and test purposes.
 * \return  number of distinct connected components
 */
size_t ExtendedCCDBG::count_connected_components(){
    std::unordered_set<unsigned> unique_set;
    for (auto &unitig : *this){
        DataAccessor<UnitigExtension>* da = unitig.getData();
        UnitigExtension* ue = da->getData(unitig);
        unique_set.insert(seqan::findSet(UF, ue->getID()));
    }
    return unique_set.size();
}


/*!
 * \fn      bool ExtendedCCDBG::connected_components(const CDBG_Build_opt &graph_options)
 * \brief   This function computes the connected components for at the current state of the graph.
 * \ref     seqan/include/seqan/misc/union_find.h
 * \return  true if successful
 */
bool ExtendedCCDBG::connected_components(const CCDBG_Build_opt &graph_options){
    // initiate UF structure
    if (graph_options.verbose) std::cout << "[VERBOSE] Initiating UNION-FIND" << std::endl;
    resize(UF, (*this).size()+1);   // however, UF needs to be +1 bigger than the graph size. Probably due to start index 1 of IDs

#ifdef DEBUG
    std::cout << "UF size " << seqan::length(UF._values) << std::endl;
#endif // DEBUG

    // run UF joining
    if (graph_options.verbose) std::cout << "[VERBOSE] Running UNION-FIND" << std::endl;
    for (auto &unitig : *this){
        // TODO: progress indicator here

        /*  Get all predecessors AND successors of a node (unitig).
         *  I need to iterate through both since either one of them could miss
         *  links, and therefore split components, where both unitigs
         *  "face each other", e.g.:
         *        u1 ----------->
         *                 <-------------- u2
         * which is in GFA:
         *      L   u1  +   u2  -
         *      L   u2  +   u1  -
         */
        DataAccessor<UnitigExtension>* da = unitig.getData();
        UnitigExtension* ue = da->getData(unitig);
        size_t unitig_id = ue->getID();

        for (auto &it_pre : unitig.getPredecessors()){
            DataAccessor<UnitigExtension>* da = it_pre.getData();
            UnitigExtension* ue = da->getData(it_pre);
            unsigned pre_id = ue->getID();
            seqan::joinSets(UF, seqan::findSet(UF, unitig_id), seqan::findSet(UF, pre_id));
        }

        for (auto &it_suc : unitig.getSuccessors()){
            DataAccessor<UnitigExtension>* da = it_suc.getData();
            UnitigExtension* ue = da->getData(it_suc);
            unsigned suc_id = ue->getID();
            seqan::joinSets(UF, seqan::findSet(UF, unitig_id), seqan::findSet(UF, suc_id));
        }
    }

    return true;
}


/*!
 * \fn          inline float ExtendedCCDBG::entropy(const std::string &sequence)
 * \brief       This function computes an entropy for a given string that can be used to filter/mark low complexity
 *              sequences. If all dimers are equaly distributed the entropy is high ("highly chaotic system"), if
 *              all dimers follow a certain pattern the entropy is low ("highly ordered system"). We'd probably like
 *              to mark low entropy unitigs since they have a chance to disrupt/branch the de Bruijn Graph.
 * \remark      Function taken from PopIns.
 * \return      The entropy [0,1] of all binucleotides.
 */
inline float ExtendedCCDBG::entropy(const std::string &sequence){
    // create a dictionary counting the occurrence of all dinucleotides
    unordered_map<std::string, unsigned> diCounts(16);
    unsigned counted = 0;
    for (unsigned i = 0; i < sequence.length()-1; ++i){
        std::string dimer = sequence.substr(i,2);
        if (sequence[i]!='N' && sequence[i+1]!='N'){
            // set if dimer not in counter table yet
            if(diCounts.find(dimer) == diCounts.end()){
                diCounts[dimer] = 1;
                counted++;
            }
            // otherwise increase
            else{
                diCounts[dimer] += 1;
                counted++;
            }
        }
    }

    // calculate the entropy for dinucleotide counts
    unsigned entropy = 0;
    for(unordered_map<std::string,unsigned>::const_iterator it = diCounts.cbegin(); it != diCounts.cend(); ++it){
        if (it->second == 0) continue;
        float p = float(it->second) / counted;
        entropy -= p * log(p) / log(2);
    }

    return entropy / 4;
}


inline void ExtendedCCDBG::DFS_cleaner(){
    for (auto &ucm : *this){
        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* ue = da->getData(ucm);

        ue->set_undiscovered_fw();
        ue->set_undiscovered_bw();
    }
}


inline void ExtendedCCDBG::DFS_cleaner_seen_only(){
    for (auto &ucm : *this){
        DataAccessor<UnitigExtension>* da = ucm.getData();
        UnitigExtension* ue = da->getData(ucm);

        // reset only internal nodes
        if (ue->is_seen_fw())
            ue->set_undiscovered_fw();
        if (ue->is_seen_bw())
            ue->set_undiscovered_bw();
    }
}


/*!
 * \fn      inline uint8_t ExtendedCCDBG::whereToGo(const UnitigColorMap< UnitigExtension >& um, const UnitigColorMap< UnitigExtension >& src) const
 * \brief   This function tests the predecessors P of a unitig u. If the searched unitig src is in P, then
 *          whereToGo() returns GO_FORWARD, denoting the traversal has to continue in the successors of u. If src is
 *          not in P, then whereToGo() returns GO_BACKWARD, denoting the traversal has to continue in the predecessors
 *          of u.
 * \details NOTE: This function is a lowest-level indication where to go, independent of the unitig's orientation. If
 *          the orientation of a unitig is rev-comp, then the result might be GO_BACKWARD while we still consider it
 *          a forward motion with respect to the traversal.
 * \return  DIRECTION
 */
inline uint8_t ExtendedCCDBG::whereToGo(const UnitigColorMap< UnitigExtension >& um, const UnitigColorMap< UnitigExtension >& src) const{
    uint8_t ret = GO_BACKWARD;
    for (auto &predecessor : um.getPredecessors())
         if (predecessor == src)    /* NOTE: If traversal is within a small-loop, it will always evaluate to true. Therefore loop cases have to be catched during traversal.*/
             ret = GO_FORWARD;
    return ret;
}


/*!
 * \fn      inline uint8_t ExtendedCCDBG::whereFrom(const UnitigColorMap< UnitigExtension >& um, const UnitigColorMap< UnitigExtension >& src) const
 * \brief   Reverses the direction of whereToGo(). E.g. if the answer is "GO_BACKWARD", then
 *          this function indicates where to go to reach the source (src).
 *          src  ----------
 *          um          --------->
 * \details This function can also be used to detemine a neighbor's (NBR) orientation with respect to the current unitig (CU) if src=CU and um=NBR.
 * \return  DIRECTION
 */
inline uint8_t ExtendedCCDBG::whereFrom(const UnitigColorMap< UnitigExtension >& um, const UnitigColorMap< UnitigExtension >& src) const{
    return (whereToGo(um, src) == GO_FORWARD) ? GO_BACKWARD : GO_FORWARD;
}


/*!
 * \fn      Traceback ExtendedCCDBG::DFS_Init(const UnitigColorMap<UnitigExtension> &ucm, Setcover<> &sc, const bool verbose)
 * \brief   This function initiates the recursion of the directed DFS.
 * \return  Traceback object
 */
Traceback ExtendedCCDBG::DFS_Init(const UnitigColorMap<UnitigExtension> &ucm,
                                  Setcover<> &sc,
                                  const bool verbose,
                                  const unsigned max_paths){
    Traceback tb;
    tb.recursion_priority_counter = max_paths - 1;

    // get data of current unitig
    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* ue = da->getData(ucm);
    if (verbose) cout << "I am starting at " << ue->getID() << "." << endl;

    // abort if start node was traversed already (can happen when max_paths > 1)
    if (ue->is_visited_fw()){     // could be fw or bw, doesn't matter
        if (verbose) cout << "I abort at " << ue->getID() << ", node was visited already." << endl;
        return tb;
    }

    BackwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> bw_neighbors = ucm.getPredecessors();
    ForwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> fw_neighbors = ucm.getSuccessors();

    // ---------------------
    // | handle singletons |
    // ---------------------
    if (!bw_neighbors.hasPredecessors() && !fw_neighbors.hasSuccessors()){

        VSequences vseqs;
        ucm.strand ? vseqs.push_back(ucm.referenceUnitigToString()) : vseqs.push_back(reverse_complement(ucm.referenceUnitigToString()));
        tb.push_back(vseqs);
        tb.recursive_return_status = true;

        sc.add(ue->getID());
        sc.unify();
        sc.clear();

        return tb;
    }

    // -------------------------
    // | handle internal nodes |
    // -------------------------
    else if(bw_neighbors.hasPredecessors() && fw_neighbors.hasSuccessors()){
        // NOTE: just do nothing in this case (internal node)
        if (verbose) cout << "Returning from internal node." << endl;
        return tb;
    }
    // -------------------
    // | start traversal |
    // -------------------
    else{
        // yShape check
        if (yStemCheck(ucm, verbose)){
            if (verbose) cout << "Discarding startnode " << ue->getID() << " since it is a heavy Y-chain element." << endl;
            return tb;
        }

        const uint8_t direction = !bw_neighbors.hasPredecessors() ? GO_FORWARD : GO_BACKWARD;
        size_t nb_colors = this->getNbColors();
        std::vector<bool> start_vec(nb_colors, false);

        // ----------------
        // | go backward  |
        // ----------------
        if (direction==GO_BACKWARD){
            if (verbose) cout << "I am setting " << ue->getID() << " to seen (bw)." << endl;
            ue->set_seen_bw();

            // init start_vec bw traversal
            const const_UnitigColorMap<UnitigExtension> ucm_head_kmer = ucm.getKmerMapping(0);
            const UnitigColors* ucm_head_uc = ucm_head_kmer.getData()->getUnitigColors(ucm_head_kmer);
            for (size_t color_id = 0; color_id != nb_colors; ++color_id)
                if (ucm_head_uc->contains(ucm_head_kmer, color_id))
                    start_vec[color_id] = true;

            // =====================================
            // | checkpoint: color ranked neighbors |
            // =====================================
            /**/
            // get pairs (rel. overlap, ID) of all neighbors in descending order
            std::multimap<float, unsigned, GreaterThan> descendingSortedNeighbors;              // IDEA: for more decision criteria, encode this in the compare functor
            sortNeighbors(bw_neighbors, start_vec, descendingSortedNeighbors);
            /**/

            // =====================================
            // | checkpoint: color ranked neighbors |
            // =====================================
            /**/
            for (auto it = descendingSortedNeighbors.cbegin(); it != descendingSortedNeighbors.cend(); ++it){
            /**/

                if (tb.recursion_priority_counter < max_paths){             // NOTE: add here " && tb.recursive_return_status==false" to restrict DFS to one traceback path only

                    for (auto &neighbor : bw_neighbors){
                        // =====================================
                        // | checkpoint: color ranked neighbors |
                        // =====================================
                        /**/
                        DataAccessor<UnitigExtension>* neighbor_da = neighbor.getData();
                        UnitigExtension* neighbor_ue = neighbor_da->getData(neighbor);
                        unsigned ue_id = neighbor_ue->getID();
                        unsigned neighbor_id = it->second;

                        if (ue_id == neighbor_id){
                        /**/
                            DFS_case(ucm, neighbor, start_vec, tb, sc, verbose, max_paths);
                        // =====================================
                        // | checkpoint: color ranked neighbors |
                        // =====================================
                        /**/
                            break;  // since only one neighbor ucm will match the n-th best neighbor ID, we can break after we found it
                        }
                        /**/

                    }
                }
                // =====================================
                // | checkpoint: color ranked neighbors |
                // =====================================
                /**/
                else{break;}
            }
                /**/
        }

        // ---------------
        // | go forward  |
        // ---------------
        else{   // direction==GO_FORWARD
            if (verbose) cout << "I am setting " << ue->getID() << " to seen (fw)." << endl;
            ue->set_seen_fw();

            // init start_vec fw traversal
            const const_UnitigColorMap<UnitigExtension> ucm_tail_kmer = ucm.getKmerMapping(ucm.size - static_cast<size_t>(this->getK()));
            const UnitigColors* ucm_tail_uc = ucm_tail_kmer.getData()->getUnitigColors(ucm_tail_kmer);
            for (size_t color_id = 0; color_id != nb_colors; ++color_id)
                if (ucm_tail_uc->contains(ucm_tail_kmer, color_id))
                    start_vec[color_id] = true;

            // =====================================
            // | checkpoint: color ranked neighbors |
            // =====================================
            /**/
            // get pairs (rel. overlap, ID) of all neighbors in descending order
            std::multimap<float, unsigned, GreaterThan> descendingSortedNeighbors;              // IDEA: for more decision criteria, encode this in the compare functor
            sortNeighbors(fw_neighbors, start_vec, descendingSortedNeighbors);
            /**/

            // =====================================
            // | checkpoint: color ranked neighbors |
            // =====================================
            /**/
            for (auto it = descendingSortedNeighbors.cbegin(); it != descendingSortedNeighbors.cend(); ++it){
            /**/

                if (tb.recursion_priority_counter < max_paths){             // NOTE: add here " && tb.recursive_return_status==false" to restrict DFS to one traceback path only

                    for (auto &neighbor : fw_neighbors){
                        // =====================================
                        // | checkpoint: color ranked neighbors |
                        // =====================================
                        /**/
                        DataAccessor<UnitigExtension>* neighbor_da = neighbor.getData();
                        UnitigExtension* neighbor_ue = neighbor_da->getData(neighbor);
                        unsigned ue_id = neighbor_ue->getID();
                        unsigned neighbor_id = it->second;

                        if (ue_id == neighbor_id){
                        /**/
                            DFS_case(ucm, neighbor, start_vec, tb, sc, verbose, max_paths);
                        // =====================================
                        // | checkpoint: color ranked neighbors |
                        // =====================================
                        /**/
                            break;  // since only one neighbor ucm will match the n-th best neighbor ID, we can break after we found it
                        }
                        /**/

                    }
                }
                // =====================================
                // | checkpoint: color ranked neighbors |
                // =====================================
                /**/
                else{break;}
            }
                /**/
        }

    }

    sc.clear();     // NOTE: sc.del() would be sufficient here I guess.

    if (tb.recursive_return_status){        // mark start nodes only if at least one traceback happened. Important for #recursion attempts (c) > 1
        ue->set_visited_fw();
        ue->set_visited_bw();
        if (verbose) cout << "I am setting " << ue->getID() << " to visited (both)." << endl;
        if (verbose) cout << "I am done with " << ue->getID() << endl;

        sc.addStartnode(ue->getID());
    }

    return tb;
}


/*!
 * \fn      Traceback ExtendedCCDBG::DFS_Visit(const UnitigColorMap<UnitigExtension> &ucm,
 *                                             std::vector<bool> &start_vec,
 *                                             const uint8_t src_direction,
 *                                             Setcover<> &sc,
 *                                             const bool verbose)
 * \brief   This function executes the recursion of the directed DFS.
 * \param   ucm is the current node during traversal
 * \param   src_direction is the direction you'd need to go if you wanted to go back from ucm to the previous node of the traversal.
 *          Therefore, for further traversal you need to follow the opposite direction.
 * \param   start_direction is the initial traversal direction from the DFS' start_ucm to the 2nd unitig of the traversal
 * \return  Traceback object
 */
Traceback ExtendedCCDBG::DFS_Visit(const UnitigColorMap<UnitigExtension> &ucm,
                                   std::vector<bool> &start_vec,
                                   const uint8_t src_direction,
                                   Setcover<> &sc,
                                   const bool verbose,
                                   const unsigned max_paths){
    Traceback tb;

    // get data of current unitig
    DataAccessor<UnitigExtension>* da = ucm.getData();
    UnitigExtension* ue = da->getData(ucm);

    uint8_t traversal_direction = (src_direction==GO_BACKWARD) ? GO_FORWARD : GO_BACKWARD;  // NOTE: I guess this could be avoided if I'd pass traversal direction directly in DFS_visit

    // keep start_vec up-to-date
    update_start_vec(start_vec, ucm);

    // --------------------------
    // |  if traverse backward  |
    // --------------------------
    if (traversal_direction==GO_BACKWARD){
        if (verbose) cout << "I am setting " << ue->getID() << " to seen (bw)." << endl;
        ue->set_seen_bw();  // NOTE: I keep this marked since we (in it's current version) only report one path per source-sink note

        // -----------------------
        // |  if stop by colors  |  // TODO: this should never be the case now, delete?
        // -----------------------
        if (is_empty_start_vec(start_vec)){
            if (verbose) cout << "I see " << ue->getID() << " does not satisfy the color criteria." << endl;
            if (verbose) cout << "Traversal will not go further here." << endl;

            return tb;
        }

        // ------------------
        // |  if sink node  |
        // ------------------
        BackwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> bw_neighbors = ucm.getPredecessors();

        if ( !bw_neighbors.hasPredecessors() && !ue->is_visited_bw() ){   // visited check could be bw/fw; visited check necessary to avoid RC path
            if (verbose) cout << "I see " << ue->getID() << " has no predecessors and is not visited." << endl;
            if (verbose) cout << "I will trigger traceback from " << ue->getID() << endl;

            VSequences vseqs;
            ucm.strand ? vseqs.push_back(ucm.referenceUnitigToString()) : vseqs.push_back(reverse_complement(ucm.referenceUnitigToString()));
            tb.push_back(vseqs);

            tb.recursive_return_status = true;
            tb.recursion_priority_counter += 1;
            if (verbose) cout << "check 1" << endl;

            sc.add(ue->getID());
            if (verbose) cout << "check 2" << endl;

            sc.print_current();
            sc.unify();
            if (verbose) cout << "check 3" << endl;

            return tb;
        }

        // -------------------------
        // |  if traverse further  |
        // -------------------------

        // =====================================
        // | checkpoint: color ranked neighbors |
        // =====================================
        /**/
        // get pairs (rel. overlap, ID) of all neighbors in descending order
        std::multimap<float, unsigned, GreaterThan> descendingSortedNeighbors;              // IDEA: for more decision criteria, encode this in the compare functor
        sortNeighbors(bw_neighbors, start_vec, descendingSortedNeighbors);
        /**/

        // =====================================
        // | checkpoint: color ranked neighbors |
        // =====================================
        /**/
        for (auto it = descendingSortedNeighbors.cbegin(); it != descendingSortedNeighbors.cend(); ++it){
        /**/

            if (tb.recursion_priority_counter < max_paths){             // NOTE: add here " && tb.recursive_return_status==false" to restrict DFS to one traceback path only

                for (auto &neighbor : bw_neighbors){
                    // =====================================
                    // | checkpoint: color ranked neighbors |
                    // =====================================
                    /**/
                    DataAccessor<UnitigExtension>* neighbor_da = neighbor.getData();
                    UnitigExtension* neighbor_ue = neighbor_da->getData(neighbor);
                    unsigned ue_id = neighbor_ue->getID();
                    unsigned neighbor_id = it->second;

                    if (ue_id == neighbor_id){
                    /**/
                        DFS_case(ucm, neighbor, start_vec, tb, sc, verbose, max_paths);
                    // =====================================
                    // | checkpoint: color ranked neighbors |
                    // =====================================
                    /**/
                        break;  // since only one neighbor ucm will match the n-th best neighbor ID, we can break after we found it
                    }
                    /**/

                }
            }
            // =====================================
            // | checkpoint: color ranked neighbors |
            // =====================================
            /**/
            else{break;}
        }
            /**/
    }

    // -------------------------
    // |  if traverse forward  |
    // -------------------------
    else{   // traversal_direction==GO_FORWARD
        if (verbose) cout << "I am setting " << ue->getID() << " to seen (fw)." << endl;
        ue->set_seen_fw();  // NOTE: I keep this marked since we (in it's current version) only report one path per source-sink note

        // -----------------------
        // |  if stop by colors  |  // TODO: this should never be the case now, delete?
        // -----------------------
        if (is_empty_start_vec(start_vec)){
            if (verbose) cout << "I see " << ue->getID() << " does not satisfy the color criteria." << endl;
            if (verbose) cout << "Traversal will not go further here." << endl;

            return tb;
        }

        // ------------------
        // |  if sink node  |
        // ------------------
        ForwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> fw_neighbors = ucm.getSuccessors();

        if ( !fw_neighbors.hasSuccessors() && !ue->is_visited_fw() ){   // visited check could be bw/fw; visited check necessary to avoid RC path
            if (verbose) cout << "I see " << ue->getID() << " has no successors and is not visited." << endl;
            if (verbose) cout << "I will trigger traceback from " << ue->getID() << endl;

            VSequences vseqs;
            ucm.strand ? vseqs.push_back(ucm.referenceUnitigToString()) : vseqs.push_back(reverse_complement(ucm.referenceUnitigToString()));
            tb.push_back(vseqs);

            tb.recursive_return_status = true;
            tb.recursion_priority_counter += 1;
            if (verbose) cout << "check 1" << endl;

            sc.add(ue->getID());
            if (verbose) cout << "check 2" << endl;

            sc.print_current();
            sc.unify();
            if (verbose) cout << "check 3" << endl;

            return tb;
        }

        // -------------------------
        // |  if traverse further  |
        // -------------------------

        // =====================================
        // | checkpoint: color ranked neighbors |
        // =====================================
        /**/
        // get pairs (rel. overlap, ID) of all neighbors in descending order
        std::multimap<float, unsigned, GreaterThan> descendingSortedNeighbors;              // IDEA: for more decision criteria, encode this in the compare functor
        sortNeighbors(fw_neighbors, start_vec, descendingSortedNeighbors);
        /**/

        // =====================================
        // | checkpoint: color ranked neighbors |
        // =====================================
        /**/
        for (auto it = descendingSortedNeighbors.cbegin(); it != descendingSortedNeighbors.cend(); ++it){
        /**/

            if (tb.recursion_priority_counter < max_paths){             // NOTE: add here " && tb.recursive_return_status==false" to restrict DFS to one traceback path only

                for (auto &neighbor : fw_neighbors){
                    // =====================================
                    // | checkpoint: color ranked neighbors |
                    // =====================================
                    /**/
                    DataAccessor<UnitigExtension>* neighbor_da = neighbor.getData();
                    UnitigExtension* neighbor_ue = neighbor_da->getData(neighbor);
                    unsigned ue_id = neighbor_ue->getID();
                    unsigned neighbor_id = it->second;

                    if (ue_id == neighbor_id){
                    /**/
                        DFS_case(ucm, neighbor, start_vec, tb, sc, verbose, max_paths);
                    // =====================================
                    // | checkpoint: color ranked neighbors |
                    // =====================================
                    /**/
                        break;  // since only one neighbor ucm will match the n-th best neighbor ID, we can break after we found it
                    }
                    /**/

                }
            }
            // =====================================
            // | checkpoint: color ranked neighbors |
            // =====================================
            /**/
            else{break;}
        }
            /**/
    }

    /* Whenever the traversal jumps back from a deeper recursion level, the current ucm id has
       to be deleted from the current path of the setcover instance, independed whether the
       traceback is happening or not. */
    sc.del();

    return tb;
}


/*!
 * \fn      void ExtendedCCDBG::DFS_case(const UnitigColorMap<UnitigExtension> &ucm,
 *                                       const UnitigColorMap<UnitigExtension> &neighbor,
 *                                       std::vector<bool> &start_vec,
 *                                       Traceback &tb,
 *                                       Setcover<> &sc,
 *                                       const bool verbose)
 * \brief   This function contains the code to determine the neighbor's (pre/suc) orientation, the recursion call of the according case and
 *          the management of a (recursively) returned traceback instance.
 * \param   start_ucm is the DFS' initial unitig
 * \param   start_direction is the initial traversal direction from the DFS' start_ucm to the 2nd unitig of the traversal
 * \return  void
 */
void ExtendedCCDBG::DFS_case(const UnitigColorMap<UnitigExtension> &ucm,
                             const UnitigColorMap<UnitigExtension> &neighbor,
                             std::vector<bool> &start_vec,
                             Traceback &tb,
                             Setcover<> &sc,
                             const bool verbose,
                             const unsigned max_paths){

    DataAccessor<UnitigExtension>* ucm_da = ucm.getData();
    UnitigExtension* ucm_ue = ucm_da->getData(ucm);

    DataAccessor<UnitigExtension>* neighbor_da = neighbor.getData();
    UnitigExtension* neighbor_ue = neighbor_da->getData(neighbor);

    // ----------------
    // |   Case 1+3   |
    // ----------------
    if (whereFrom(neighbor, ucm)==GO_BACKWARD){
        /*  Case 3:                             Case 1:
        *   SRC ------->                OR              -------> SRC
        *           -------> SUC                PRE <--------
        */
        if (neighbor_ue->is_undiscovered_fw()){
            // ---------------
            // |  Recursion  |
            // ---------------

            if (verbose) cout << "I am at " << ucm_ue->getID() << " and will go forward to " << neighbor_ue->getID() << endl;

            sc.add(ucm_ue->getID());

            Traceback returned_tb = DFS_Visit(neighbor, start_vec, GO_BACKWARD, sc, verbose, max_paths);

            if (verbose) cout << "I jumped back to ID " << ucm_ue->getID() << endl;

            // -----------------------------------------
            // |  check recursively returned tb object |
            // -----------------------------------------
            if (returned_tb.recursive_return_status){

                for (auto it = returned_tb.begin(); it != returned_tb.end(); ++it)
                    ucm.strand ? it->push_back(ucm.referenceUnitigToString()) : it->push_back(reverse_complement(ucm.referenceUnitigToString()));
                tb.recursive_return_status = true;
                tb.join(returned_tb);

            }
        }
        else{
            // This means I have see nthis unitig before. Do nothing at the moment. Happens e.g. at loops.
            if (verbose) cout << "I have seen " << neighbor_ue->getID() << " before. I'll go back without traceback." << endl;
        }
    }
    // ----------------
    // |   Case 2+4   |
    // ----------------
    else{   // whereFrom(neighbor, ucm)==GO_FORWARD
        /*  Case 4:                             Case 2:
        *   SRC ------->                OR              -------> SRC
        *           <------- SUC                PRE ------->
        */
        if (neighbor_ue->is_undiscovered_bw()){
            // ---------------
            // |  Recursion  |
            // ---------------

            if (verbose) cout << "I am at " << ucm_ue->getID() << " and will go forward to " << neighbor_ue->getID() << endl;

            sc.add(ucm_ue->getID());

            Traceback returned_tb = DFS_Visit(neighbor, start_vec, GO_FORWARD, sc, verbose, max_paths);

            if (verbose) cout << "I jumped back to ID " << ucm_ue->getID() << endl;

            // -----------------------------------------
            // |  check recursively returned tb object |
            // -----------------------------------------
            if (returned_tb.recursive_return_status){

                for (auto it = returned_tb.begin(); it != returned_tb.end(); ++it)
                    ucm.strand ? it->push_back(ucm.referenceUnitigToString()) : it->push_back(reverse_complement(ucm.referenceUnitigToString()));
                tb.recursive_return_status = true;
                tb.join(returned_tb);

            }
        }
        else{
            // This means I have see nthis unitig before. Do nothing at the moment. Happens e.g. at lops.
            if (verbose) cout << "I have seen " << neighbor_ue->getID() << " before. I'll go back without traceback." << endl;
        }
    }

}


/*!
 * \fn      inline bool ExtendedCCDBG::endsHaveCommonColor(const UnitigColorMap<UnitigExtension> &start_ucm,
 *                                                         const UnitigColorMap<UnitigExtension> &ucm,
 *                                                         const uint8_t start_direction,
 *                                                         const bool verbose) const
 * \brief   This function intersects the last color vector (depending on orientation) of the DFS start unitig (start_ucm) with the intersection of
 *          the current unitig's (ucm) start and end color vectors. If the outer intersection is not empty, this function returns true.
 * \details ColVec_start_? INTERSECT (ColVec_current_head INTERSECT ColVec_current_tail)
 * \return  bool, true if start and current unitig have common sample (color)
 */
inline bool ExtendedCCDBG::haveCommonColor(const UnitigColorMap<UnitigExtension> &start_ucm,
                                           const UnitigColorMap<UnitigExtension> &ucm,
                                           const uint8_t start_direction,
                                           const bool verbose) const{
    bool rValue = false;
    size_t nb_colors = this->getNbColors();

    // current unitig color vectors
    const const_UnitigColorMap<UnitigExtension> ucm_head_kmer = ucm.getKmerMapping(0);
    const const_UnitigColorMap<UnitigExtension> ucm_tail_kmer = ucm.getKmerMapping(ucm.size - static_cast<size_t>(this->getK()));
    const UnitigColors* ucm_head_uc = ucm_head_kmer.getData()->getUnitigColors(ucm_head_kmer);
    const UnitigColors* ucm_tail_uc = ucm_tail_kmer.getData()->getUnitigColors(ucm_tail_kmer);

    // DFS start node color vector
    const const_UnitigColorMap<UnitigExtension> start_kmer = (start_direction==GO_FORWARD) ? start_ucm.getKmerMapping(start_ucm.size - static_cast<size_t>(this->getK())) : start_ucm.getKmerMapping(0);
    const UnitigColors* start_uc = start_kmer.getData()->getUnitigColors(start_kmer);

    if (verbose){
        DataAccessor<UnitigExtension>* da_start_ucm = start_ucm.getData();
        UnitigExtension* ue_start_ucm = da_start_ucm->getData(start_ucm);

        DataAccessor<UnitigExtension>* da_ucm = ucm.getData();
        UnitigExtension* ue_ucm = da_ucm->getData(ucm);

        cout << "Color Compare: Current Unitig: " << ue_ucm->getID() << ", Start Unitig: " << ue_start_ucm->getID() << endl;
    }

    if (verbose)
        cout << "u_s: | u_e: | s__: | return" << endl;

    for (size_t color_id = 0; color_id != nb_colors; ++color_id){
        bool hasColor_ucm_head = ucm_head_uc->contains(ucm_head_kmer, color_id);
        bool hasColor_ucm_tail = ucm_tail_uc->contains(ucm_tail_kmer, color_id);
        bool hasColor_start    =    start_uc->contains(start_kmer,    color_id);

        rValue = true && hasColor_ucm_head && hasColor_ucm_tail && hasColor_start;

        if (verbose)
            cout << hasColor_ucm_head << "    | " << hasColor_ucm_tail << "    | " << hasColor_start << "    | " << rValue << endl;

        if (rValue) break;
    }

    if (verbose)
        if (!rValue) cout << "Color-criterion stopped traversal." << endl;

    return rValue;
}


/*!
 * \fn      bool ExtendedCCDBG::merge(const CCDBG_Build_opt &opt)
 * \return  bool; true if successful
 */
bool ExtendedCCDBG::merge(const CCDBG_Build_opt &opt, const unsigned max_paths){

    /* sanity check */
    if (!this->is_id_init())
        return false;

    std::string sv_filename = "contigs.fa";
    ofstream ofs(sv_filename, std::ofstream::out);
    if (!ofs.is_open()){
        cerr << "Error: Couldn't open ofstream for contig file." << endl;
        return false;
    }

    // =======================================
    // | checkpoint: color ranked start nodes |
    // =======================================
    /*
    std::multimap<unsigned, unsigned, GreaterThan> start_nodes;     // key: #colors, value: unitig-id
    getSourceNodes(start_nodes);
    for (auto it = start_nodes.cbegin(); it != start_nodes.cend(); ++it)
        std::cout << " [" << it->first << ':' << it->second << ']';
    std::cout << '\n';
    */

    Setcover<> sc;

    size_t sv_counter = 0;

    // max attempts per start node
    for (unsigned c = max_paths; c >= 1; --c){
        cout << "max_paths: " << c << endl;

        // =======================================
        // | checkpoint: color ranked start nodes |
        // =======================================
        /*
        for (auto it = start_nodes.cbegin(); it != start_nodes.cend(); ++it){
        */

            for (auto &unitig : *this){
                // =======================================
                // | checkpoint: color ranked start nodes |
                // =======================================
                /*
                DataAccessor<UnitigExtension>* unitig_da = unitig.getData();
                UnitigExtension* unitig_ue = unitig_da->getData(unitig);
                unsigned id = unitig_ue->getID();
                if(it->second == id){
                    std::cout << "Start with id [" << it->second << ']' << std::endl;
                */
                    if (opt.verbose) cout << " -------------------------------- " << endl;
                    Traceback tb = DFS_Init(unitig, sc, opt.verbose, max_paths);
                    if (tb.recursive_return_status){

                        if (!tb.write(ofs, opt.k, sv_counter))
                            return false;
                    }
                    DFS_cleaner_seen_only();
                //}
            }
        //}
    }

    ofs.close();

    sc.print_CSV();

    return true;
}


/*!
 * \fn      void ExtendedCCDBG::getSourceNodes(TContainer &m) const
 * \details This function fills the given container with start node IDs ranked by their
 *          amount of colors, in descending order.
 */
template <class TContainer>
void ExtendedCCDBG::getSourceNodes(TContainer &m) const {
    size_t nb_colors = this->getNbColors();

    for (auto &ucm : *this){
        BackwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, true> bw_neighbors = ucm.getPredecessors();
        ForwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, true> fw_neighbors = ucm.getSuccessors();

        bool pre = bw_neighbors.hasPredecessors();
        bool suc = fw_neighbors.hasSuccessors();

        if ((pre && !suc) || (!pre && suc)){
            // get ID
            const DataAccessor<UnitigExtension>* da = ucm.getData();
            const UnitigExtension* ue = da->getData(ucm);
            unsigned id = ue->getID();

            // get #colors in traversal direction
            const const_UnitigColorMap<UnitigExtension> start_kmer = (!pre && suc) ? ucm.getKmerMapping(ucm.size - static_cast<size_t>(this->getK())) : ucm.getKmerMapping(0);
            const UnitigColors* start_uc = start_kmer.getData()->getUnitigColors(start_kmer);
            size_t nb_colors_in_start_kmer = 0;
            for (size_t color_id = 0; color_id != nb_colors; ++color_id)
                if (start_uc->contains(start_kmer, color_id)) ++nb_colors_in_start_kmer;

            // save (#colors, ID) in multimap
            m.insert(std::pair<size_t, size_t>(nb_colors_in_start_kmer, id));
        }
    }
}


/*!
 * \fn      inline float ExtendedCCDBG::equalColorbitsRate(const std::vector<bool> &v,
 *                                                         const UnitigColorMap<UnitigExtension> &neighbor) const
 * \details This function computes the rate of equal color bit, i.e. how many times we observe
 *          0&&0 or 1&&1 in the comparison between a neighbor's colors (neighbor) and a reference
 *          vector (v) divided by the #samples in the graph.
 *          This value is used to rank neighbors for further traversal.
 * \return  float
 */
inline float ExtendedCCDBG::equalColorbitsRate(const std::vector<bool> &v,
                                               const UnitigColorMap<UnitigExtension> &neighbor) const{
    unsigned equalColorbits = 0;
    size_t nb_colors = this->getNbColors();

    // current unitig color vectors
    const const_UnitigColorMap<UnitigExtension> ucm_head_kmer = neighbor.getKmerMapping(0);
    const const_UnitigColorMap<UnitigExtension> ucm_tail_kmer = neighbor.getKmerMapping(neighbor.size - static_cast<size_t>(this->getK()));
    const UnitigColors* ucm_head_uc = ucm_head_kmer.getData()->getUnitigColors(ucm_head_kmer);
    const UnitigColors* ucm_tail_uc = ucm_tail_kmer.getData()->getUnitigColors(ucm_tail_kmer);

    for (size_t color_id = 0; color_id != nb_colors; ++color_id){
        bool hasColor_ucm_head = ucm_head_uc->contains(ucm_head_kmer, color_id);
        bool hasColor_ucm_tail = ucm_tail_uc->contains(ucm_tail_kmer, color_id);

        if (hasColor_ucm_head==hasColor_ucm_tail && hasColor_ucm_tail==v[color_id])
            ++equalColorbits;
    }

    return static_cast<float>(equalColorbits)/nb_colors;
}


/*!
 * \fn      inline void ExtendedCCDBG::update_start_vec(std::vector<bool> &start_vec,
                                                        const UnitigColorMap<UnitigExtension> &ucm) const
 * \details This function updates the color bits in a referenced vector (start_vec) with
 *          the color bits of the head and tail kmer of ucm. No return value needed.
 */
inline void ExtendedCCDBG::update_start_vec(std::vector<bool> &start_vec,
                                            const UnitigColorMap<UnitigExtension> &ucm) const {
    size_t nb_colors = this->getNbColors();

    // ucm color vectors
    const const_UnitigColorMap<UnitigExtension> ucm_head_kmer = ucm.getKmerMapping(0);
    const const_UnitigColorMap<UnitigExtension> ucm_tail_kmer = ucm.getKmerMapping(ucm.size - static_cast<size_t>(this->getK()));
    const UnitigColors* ucm_head_uc = ucm_head_kmer.getData()->getUnitigColors(ucm_head_kmer);
    const UnitigColors* ucm_tail_uc = ucm_tail_kmer.getData()->getUnitigColors(ucm_tail_kmer);

    for (size_t color_id = 0; color_id != nb_colors; ++color_id){
        bool hasColor_ucm_head = ucm_head_uc->contains(ucm_head_kmer, color_id);
        bool hasColor_ucm_tail = ucm_tail_uc->contains(ucm_tail_kmer, color_id);

        start_vec[color_id] = start_vec[color_id] &&
                              hasColor_ucm_head &&
                              hasColor_ucm_tail;
    }
}


/*!
 * \fn      inline bool ExtendedCCDBG::is_empty_start_vec(const std::vector<bool> &start_vec) const
 * \details This function checks, if a vector of bool contains any set bits (1s).
 * \return  bool; true if no entry is true
 */
inline bool ExtendedCCDBG::is_empty_start_vec(const std::vector<bool> &start_vec) const{
    if (std::find(start_vec.cbegin(), start_vec.cend(), true) == start_vec.cend())
        return true;
    return false;
}


/*!
 * \fn      yStemCheck()
 * \details This function checks whether a start node is the stem of Y-shaped connected component (CC).
 *          Speaking in terms of an antibody, the traversal of a y-shaped CC should not start in the heavy domain. If it
 *          does, there is no chance for the second (not yet traversed) light chain to end up in a traversed path since the
 *          heavy domain will be marked as visited already. To avoid this, we simply force the traversal not to start in
 *          the heavy domain.
 * \return  bool; true if startnode is Y-stem
 */
inline bool ExtendedCCDBG::yStemCheck(const UnitigColorMap<UnitigExtension> &startnode, const bool verbose){
    // TODO: can be optimized to not calculate the neighbors twice (done in DFS_Init already)
    BackwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> bw_neighbors = startnode.getPredecessors();
    ForwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> fw_neighbors = startnode.getSuccessors();

    bool isYstem = true;
    if (bw_neighbors.hasPredecessors()){
        for (auto &neighbors : bw_neighbors){
            BackwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> pre_pre = neighbors.getPredecessors();
            ForwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> pre_suc = neighbors.getSuccessors();
            size_t c_pre_pre = pre_pre.cardinality();
            size_t c_pre_suc = pre_suc.cardinality();
            if (c_pre_pre + c_pre_suc != 1){
                isYstem = false;
                break;
            }
        }
    }
    else{
        for (auto &neighbors : fw_neighbors){
            BackwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> suc_pre = neighbors.getPredecessors();
            ForwardCDBG<DataAccessor<UnitigExtension>, DataStorage<UnitigExtension>, false> suc_suc = neighbors.getSuccessors();
            size_t c_suc_pre = suc_pre.cardinality();
            size_t c_suc_suc = suc_suc.cardinality();
            if (c_suc_pre + c_suc_suc != 1){
                isYstem = false;
                break;
            }
        }
    }
    return isYstem;
}


/*!
 * \fn      sortNeighbors()
 * \details Sort the neighbors given via the NeighborCDBG and store them in the new order in TMap (IDs only, cannot save unitigs)
 * \param   TNeighborCDBG is either a ForwardCDBG or BackwardCDBG
 * \param   TMap is a sorted map with a float as key and unitig ID as value
 */
template <class TNeighborCDBG, class TMap>
inline void ExtendedCCDBG::sortNeighbors(const TNeighborCDBG &neighbors, const std::vector<bool> &start_vec, TMap &container) const{
    for (auto &neighbor : neighbors){
        DataAccessor<UnitigExtension>* neighbor_da = neighbor.getData();
        UnitigExtension* neighbor_ue = neighbor_da->getData(neighbor);
        unsigned id = neighbor_ue->getID();
        float ecr = equalColorbitsRate(start_vec, neighbor);
        if (ecr > 0.0f)
            container.insert(std::pair<float, unsigned>(ecr, id));
    }
}










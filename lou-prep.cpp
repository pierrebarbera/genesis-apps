/*
    Copyright (C) 2018 Pierre Barbera

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Pierre Barbera <pierre.barbera@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "genesis/genesis.hpp"

#include <string>

using namespace genesis;
using namespace genesis::sequence;
using namespace genesis::tree;
using namespace genesis::utils;

void fail( bool const fail_if_true, std::string const& msg ) {
    if( fail_if_true ) {
        throw std::runtime_error{ msg };
    }
}

int main( int argc, char** argv )
{
    // Check if the command line contains the right number of arguments.
    if (argc != 4) {
        throw std::runtime_error(
            std::string( "Usage: " ) + argv[0] + " <newick> <fasta> <out_dir>"
        );
    }

    Options::get().allow_file_overwriting(true);
    Logging::log_to_stdout();

    std::string tree_file(argv[1]);
    std::string msa_file(argv[2]);
    std::string out_dir(argv[3]);

    out_dir = dir_normalize_path( out_dir );

    auto msa = FastaReader().read( from_file( msa_file ) );
    auto original_tree = CommonTreeNewickReader().read( from_file( tree_file ) );

    fail( is_rooted( original_tree ),
        "Tree must be unrooted" );
    fail( msa.size() != leaf_node_count( original_tree ),
        "Tree and MSA must have same number of taxa" );
    fail( not path_exists( out_dir ),
        "No such directory: " + out_dir );

    // for every leaf ID
    auto leaf_ids = leaf_node_indices( original_tree );
    for (auto id : leaf_ids) {
        // make a copy to manipulate (relatively costly but much safer than trying to undo changes)
        CommonTree tree( original_tree );

        auto& leaf_node = tree.node_at( id );
        auto& base_node = leaf_node.link().outer().node();
        assert( &leaf_node.link().next() == &leaf_node.link() );
        assert( is_leaf(leaf_node) );
        // auto& prox_node = base_node.link().outer().node();

        auto const leaf_name = leaf_node.data<CommonNodeData>().name;
        // track pendant and proximal length
        auto const pendant_length   = leaf_node.primary_edge().data<CommonEdgeData>().branch_length;
        auto const proximal_length  = base_node.primary_edge().data<CommonEdgeData>().branch_length;

        // remove leaf from tree
        delete_leaf_node( tree, leaf_node );

        // track edge that will become the leftover edge in the pruned tree
        CommonTreeEdge* attach_edge_ptr = &base_node.primary_edge();

        // remove
        delete_linear_node( tree, base_node, []( TreeEdge& r_edge, TreeEdge& d_edge ){
            r_edge.data<CommonEdgeData>().branch_length += d_edge.data<CommonEdgeData>().branch_length;
        } );

        // determine edge_num of edge where we pruned
        size_t edge_num = 0;
        for( auto const& it : postorder( tree ) ) {
            if( it.is_last_iteration() ) { continue; }

            if( &it.edge() == attach_edge_ptr ) {
                break;
            } else {
                edge_num++;
            }
        }
        fail( edge_num >= edge_count( tree ), std::to_string(edge_num) + " vs. " + std::to_string( edge_count( tree ) ) );
        assert( edge_num < edge_count( tree ) );

        // make an output dir for this leaf
        auto cur_out_dir = dir_normalize_path( out_dir + std::to_string( id ) );
        dir_create( cur_out_dir );

        // write leaf name, edge_num, pendant and distal lengths
        std::ofstream locfile( cur_out_dir + "original_location.csv" );
        locfile << "name,edge_num,pendant,proximal\n";
        locfile << leaf_name + ",";
        locfile << std::to_string( edge_num ) + ",";
        locfile << std::to_string( pendant_length ) + ",";
        locfile << std::to_string( proximal_length ) + "\n";

        // write the tree
        CommonTreeNewickWriter().to_file( tree, cur_out_dir + "pruned.newick" );

        // write the sequence, unaligned, into a phylip file
        SequenceSet seq_set;
        auto& seq = seq_set.add( *find_sequence( msa, leaf_name ) );
        remove_all_gaps( seq );

        PhylipWriter().to_file( seq_set, cur_out_dir + "query.phylip" );

        // write out ref msa
        SequenceSet ref_set( msa );
        for (auto it = ref_set.begin(); it != ref_set.end(); ++it) {
            if( it->label() == leaf_name ) {
                ref_set.remove( it );
                break;
            }
        }

        PhylipWriter().to_file( ref_set, cur_out_dir + "reference.phylip" );
    }

    return 0;
}

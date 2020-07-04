/*
    Copyright (C) 2020 Pierre Barbera

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

#include <algorithm>
#include <fstream>
#include <string>

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

void fail( bool const failure_condition, std::string const& msg )
{
  if( failure_condition ) {
    throw std::runtime_error{ msg };
  }
}

void prune_and_move_placements( Sample& sample,
                                std::string const& prune_label ) {
  auto& tree = sample.tree();
  // search for the leaf node of the given label
  auto leaf_node  = find_node( tree, prune_label );
  fail( not leaf_node, "Given label was not found in the tree." );

  auto& edge      = leaf_node->primary_edge();
  auto& base_node = leaf_node->link().outer().node();

  // find all placements that belong to this edge
  auto cur_placements = placements_per_edge( sample, edge );

  // remove leaf from tree
  delete_leaf_node( tree, *leaf_node );

  TreeEdge* new_edge;
  // remove inner node
  delete_linear_node( tree, base_node, [&new_edge]( TreeEdge& r_edge, TreeEdge& d_edge ) {
    r_edge.data< CommonEdgeData >().branch_length += d_edge.data< CommonEdgeData >().branch_length;
    new_edge = &r_edge;
  } );

  // transfer the placements to the new edge
  for( auto& placement : cur_placements ) {
    const_cast< PqueryPlacement* >(placement)->reset_edge( *new_edge );
  }
}

/**
 *  Compares the distances of placements between two files, 
 *  assuming one has a pruned underlying tree compared to the other, 
 *  where re-pruning can equalize their topologies
 */
int main( int argc, char** argv )
{
  Options::get().allow_file_overwriting( true );
  // Activate logging.
  Logging::log_to_stdout();
  Logging::details.date = true;
  Logging::details.time = true;


  // Check if the command line contains the right number of arguments.
  if( argc != 3 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + " <jplace-file> <jplace-file-pruned>" );
  }

  LOG_INFO << "Started";
  
  // In out dirs.
  auto reader = JplaceReader();

  auto lhs = reader.read( from_file( argv[ 1 ] ) );
  auto rhs = reader.read( from_file( argv[ 2 ] ) );

  // only keep pqueries that appear in both files
  filter_pqueries_intersecting_names( lhs, rhs );

  fail( lhs.empty() or rhs.empty(), "query intersection of samples was empty" );

  // filter out all but the best hits from the pqueries (simplification for our needs)
  filter_n_max_weight_placements( lhs, 1 );
  filter_n_max_weight_placements( rhs, 1 );

  // determine which taxa were pruned out
  auto lhs_taxa = node_names( lhs.tree(), true );
  auto rhs_taxa = node_names( rhs.tree(), true );

  std::vector< std::string > missing_in_rhs;
  for( auto const& lhs_tax : lhs_taxa ) {
    if( std::find(rhs_taxa.begin(), rhs_taxa.end(), lhs_tax) == rhs_taxa.end() ) {
      missing_in_rhs.push_back( lhs_tax );
    }
  }

  std::vector< std::string > missing_in_lhs;
  for( auto const& rhs_tax : rhs_taxa ) {
    if( std::find(lhs_taxa.begin(), lhs_taxa.end(), rhs_tax) == lhs_taxa.end() ) {
      missing_in_lhs.push_back( rhs_tax );
    }
  }

  fail( missing_in_lhs.empty() and missing_in_rhs.empty(), "both empty" );

  fail( (not missing_in_lhs.empty()) and (not missing_in_rhs.empty()), "both missing taxa" );

  auto& big_sample  = missing_in_lhs.empty() ? lhs : rhs;
  auto& small_sample= missing_in_lhs.empty() ? rhs : lhs;
  auto& big_tree    = big_sample.tree();
  auto& small_tree  = small_sample.tree();
  auto& to_prune    = missing_in_lhs.empty() ? missing_in_rhs : missing_in_lhs;

  // ensure we won't prune too much
  fail( to_prune.size() > (leaf_node_count( big_tree ) - 4),
        "cannot prune this many leaves, would result in less than 4 taxa." );
  
  // prune out the extra taxa from the big tree
  for( auto const& prune_label : to_prune ) {
    prune_and_move_placements( big_sample, prune_label );
  }
  // ensure we re- calculate the indices so we can do the copy later
  reset_edge_nums( big_tree );

  fail( not compatible_trees( big_tree, small_tree ),
        "Trees are not compatible after prune" );

  // postfix the big_sample pq's to disambiguate later
  std::string postfix = "_cpy";
  for( auto& pq : big_sample ) {
    pq.name_at( 0 ).name += postfix;
  }

  // copy over pqueries from one to the other (big, which is now pruned, to small)
  copy_pqueries( big_sample, small_sample );
  
  // output node distance between all pairs of duplicates
  auto node_path_lengths = node_path_length_matrix( small_sample.tree() );
  std::vector< Pquery* > done;
  for( auto const& pq : small_sample ) {
    // check if current pq was already covered from the other direction
    if( std::find( done.begin(), done.end(), &pq ) != done.end() ) {
      continue;
    }

    auto& pq_name       = pq.name_at( 0 ).name;
    auto other_pq_name  = pq_name + postfix;

    // search for the other pquery of the same name
    auto other_pq_ptr = find_pquery( small_sample, other_pq_name );
    fail( other_pq_ptr == nullptr, "Could not find other pquery" );
    auto& other_pq = *other_pq_ptr;

    LOG_INFO << pq_name << " <-> " << other_pq_name;

    auto& first_place   = pq.placement_at( 0 );
    auto& second_place  = other_pq.placement_at( 0 );


    LOG_INFO << " path distance: " << placement_path_length_distance( first_place,
                                                                      second_place,
                                                                      node_path_lengths );


    done.push_back( other_pq_ptr );
  }

  LOG_INFO << "Finished";
  return 0;
}

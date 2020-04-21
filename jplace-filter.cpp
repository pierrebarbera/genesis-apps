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

#include <algorithm>
#include <functional>
#include <iomanip>
#include <limits>
#include <tuple>
#include <utility>

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

using weight_function = std::function< double( int ) >;

std::pair< double, double > weighted_sum_bls( Subtree const& sub, size_t const radius, weight_function& weightof )
{
  double sum          = 0.0;
  double total_weight = 0.0;

  for( auto const& it : levelorder( sub ) ) {
    assert( radius <= std::numeric_limits< int >::max() );

    if( it.depth() > static_cast< int >( radius ) ) {
      break;
    }
    auto const& edge = it.edge();

    auto bl       = edge.data< CommonEdgeData >().branch_length;
    double weight = weightof( it.depth() + 1 );

    sum += weight * bl;
    total_weight += weight;
  }

  return std::make_pair( sum, total_weight );
}

double weighted_avg_bl( TreeEdge const& edge, size_t const radius, weight_function& weightof )
{
  auto bl       = edge.data< CommonEdgeData >().branch_length;
  double weight = weightof( 0 );

  double sum          = weight * bl;
  double total_weight = weight;

  // make two subtrees at each end of the starting edge
  Subtree lhs( edge.primary_node() );
  Subtree rhs( edge.secondary_node() );

  if( radius > 0 ) {
    auto lhs_result = weighted_sum_bls( lhs, radius - 1, weightof );
    auto rhs_result = weighted_sum_bls( rhs, radius - 1, weightof );

    sum += lhs_result.first + rhs_result.first;
    total_weight += lhs_result.second + rhs_result.second;
  }

  return sum / total_weight;
}

std::vector< double > get_local_bl_avg_weighted( PlacementTree const& tree, size_t const radius, weight_function& weightof )
{
  std::vector< double > bl_map( tree.edge_count() );
  for( auto const& edge : tree.edges() ) {
    auto const edge_num = edge.data< PlacementEdgeData >().edge_num();

    bl_map[ edge_num ] = weighted_avg_bl( edge, radius, weightof );
  }

  return bl_map;
}

double subtree_max_bl( Subtree const& sub, size_t const radius )
{
  double max = 0.0;

  for( auto const& it : levelorder( sub ) ) {
    assert( radius <= std::numeric_limits< int >::max() );

    if( it.depth() > static_cast< int >( radius ) ) {
      break;
    }
    auto const& edge = it.edge();

    auto bl = edge.data< CommonEdgeData >().branch_length;
    max     = std::max( max, bl );
  }

  return max;
}

double local_max_bl( TreeEdge const& edge, size_t const radius )
{
  double max = edge.data< CommonEdgeData >().branch_length;

  // make two subtrees at each end of the starting edge
  Subtree lhs( edge.primary_node() );
  Subtree rhs( edge.secondary_node() );

  if( radius > 0 ) {
    auto lhs_max = subtree_max_bl( lhs, radius - 1 );
    auto rhs_max = subtree_max_bl( rhs, radius - 1 );

    max = std::max( max, std::max( lhs_max, rhs_max ) );
  }

  return max;
}

std::vector< double > get_local_max_bl( PlacementTree const& tree, size_t const radius )
{
  std::vector< double > bl_map( tree.edge_count() );
  for( auto const& edge : tree.edges() ) {
    auto const edge_num = edge.data< PlacementEdgeData >().edge_num();

    bl_map[ edge_num ] = local_max_bl( edge, radius );
  }

  return bl_map;
}

/**
 *
 */
int main( int argc, char** argv )
{
  // Check if the command line contains the right number of arguments.
  if( argc < 2 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + " <jplace-files...>\n" );
  }

  int radius = 10; //std::stoi(argv[1]);
  radius     = ( radius < 0 ) ? std::numeric_limits< int >::max() : radius;

  double thresh_mult = 4; //std::stod(argv[2]);

  std::vector< std::string > jplace_files;
  for( int i = 1; i < argc; ++i ) {
    jplace_files.emplace_back( argv[ i ] );
  }

  size_t total_trash = 0;

  weight_function identity    = []( int ) -> double { return 1; };
  weight_function sigmoid     = []( int x ) -> double { return 1.0 / ( 1 + exp( x - 4 ) ); };
  weight_function oneoverx    = []( int x ) -> double { return 1.0 / ( x + 1 ); };
  weight_function oneoverxsqa = []( int x ) -> double { return 1.0 / pow( x + 1, 2 ); };
  weight_function oneoverxtox = []( int x ) -> double { return 1.0 / pow( x + 1, x ); };
  weight_function exponential = []( int x ) -> double { return exp( -x / 5 ); };

  for( auto& filename : jplace_files ) {
    size_t trash_count = 0;

    std::cout << "Filtering " << filename << std::endl;
    Sample sample = JplaceReader().read( from_file( filename ) );
    Sample filtered_sample( sample.tree() );

    // precompute the average local branch length
    auto bl_map_avg          = get_local_bl_avg_weighted( sample.tree(), radius, identity );
    auto bl_map_weighted_avg = get_local_bl_avg_weighted( sample.tree(), radius, oneoverx );
    auto bl_map_max          = get_local_max_bl( sample.tree(), radius );

    for( auto& pq : sample ) {
      sort_placements_by_weight( pq );
      auto& p = pq.placement_at( 0 );

      // capture blame info
      auto const edge_num = p.edge_num();

      auto name = p.edge().primary_node().data< CommonNodeData >().name;
      if( name.empty() ) {
        name = p.edge().secondary_node().data< CommonNodeData >().name;
      }

      // ==================
      // local thresholding
      // ==================
      // if (p.pendant_length > thresh_mult * bl_map_avg[edge_num]) {
      //     ++discarded_overavg;
      // }

      // if (p.pendant_length > thresh_mult * bl_map_weighted_avg[edge_num]) {
      //     ++discarded_overavg_weighted;
      // }

      if( p.pendant_length > thresh_mult * bl_map_max[ edge_num ] ) {
        ++trash_count;
      } else {
        filtered_sample.add( pq );
      }
    }

    auto filtered_name = "filtered_" + filename;

    JplaceWriter().to_file( filtered_sample, filtered_name );

    std::cout << "Done! " << trash_count << " queries removed. Output: " << filtered_name << std::endl;

    total_trash += trash_count;
  }

  std::cout << "All done! total removed: " << total_trash << std::endl;

  return 0;
}

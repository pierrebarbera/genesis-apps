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

  auto sampleset = JplaceReader().read( from_files( jplace_files ) );

  if( sampleset.empty() ) {
    throw std::invalid_argument{ "Must supply at least one valid jplace file!" };
  }

  size_t total_entries = 0;
  for( auto& sample : sampleset ) {
    total_entries += sample.size();
  }

  Matrix< double > mat( total_entries, 2 );

  size_t row = 0;

  size_t discarded_overmax          = 0;
  size_t discarded_overavg          = 0;
  size_t discarded_overavg_weighted = 0;

  weight_function identity    = []( int ) -> double { return 1; };
  weight_function sigmoid     = []( int x ) -> double { return 1.0 / ( 1 + exp( x - 4 ) ); };
  weight_function oneoverx    = []( int x ) -> double { return 1.0 / ( x + 1 ); };
  weight_function oneoverxsqa = []( int x ) -> double { return 1.0 / pow( x + 1, 2 ); };
  weight_function oneoverxtox = []( int x ) -> double { return 1.0 / pow( x + 1, x ); };
  weight_function exponential = []( int x ) -> double { return exp( -x / 5 ); };

  std::vector< double > pendant_lengths;
  pendant_lengths.resize( total_entries );

  struct blamestruct {
    size_t edge_num  = 0;
    std::string name = "";

    size_t for_overavg          = 0;
    size_t for_overavg_weighted = 0;
    size_t for_overmax          = 0;
    size_t for_stddev_5         = 0;
  };

  size_t const total_edges = sampleset[ 0 ].tree().edge_count();

  std::vector< blamestruct > blame( total_edges );

  for( auto& sample : sampleset ) {

    // precompute the average local branch length
    auto bl_map_avg          = get_local_bl_avg_weighted( sample.tree(), radius, identity );
    auto bl_map_weighted_avg = get_local_bl_avg_weighted( sample.tree(), radius, oneoverx );
    auto bl_map_max          = get_local_max_bl( sample.tree(), radius );

    for( auto& pq : sample ) {
      sort_placements_by_weight( pq );
      auto& p = pq.placement_at( 0 );

      // capture the pendant lengths
      pendant_lengths.push_back( p.pendant_length );

      // capture blame info
      auto const edge_num        = p.edge_num();
      blame[ edge_num ].edge_num = edge_num;

      auto name = p.edge().primary_node().data< CommonNodeData >().name;
      if( name.empty() ) {
        name = p.edge().secondary_node().data< CommonNodeData >().name;
      }

      if( not name.empty() ) {
        blame[ edge_num ].name = name;
      }

      // mat.at(row, 0) = p.like_weight_ratio;
      // mat.at(row, 1) = p.pendant_length;

      // ==================
      // local thresholding
      // ==================
      if( p.pendant_length > thresh_mult * bl_map_avg[ edge_num ] ) {
        ++discarded_overavg;
        blame[ edge_num ].for_overavg++;
      }

      if( p.pendant_length > thresh_mult * bl_map_weighted_avg[ edge_num ] ) {
        ++discarded_overavg_weighted;
        blame[ edge_num ].for_overavg_weighted++;
      }

      if( p.pendant_length > thresh_mult * bl_map_max[ edge_num ] ) {
        ++discarded_overmax;
        blame[ edge_num ].for_overmax++;
      }

      row++;
    }
  }

  // sort the pendant length set to enable all possible statistical tools
  std::sort( std::begin( pendant_lengths ), std::end( pendant_lengths ) );

  std::cout << "~~~ Program Settings ~~~\n";
  std::cout << "Locality radius:\t" << radius << "\n";
  std::cout << "Threshold multiplier:\t" << thresh_mult << "\n";

  std::cout << "\n~~~ Basic Info ~~~\n";
  std::cout << "Number of input files:\t" << sampleset.size() << "\n";
  std::cout << "Number of queries:\t" << total_entries << "\n";

  auto minmax     = finite_minimum_maximum( std::begin( pendant_lengths ), std::end( pendant_lengths ) );
  auto meanstddev = mean_stddev( std::begin( pendant_lengths ), std::end( pendant_lengths ) );
  auto med        = median( std::begin( pendant_lengths ), std::end( pendant_lengths ) );

  std::cout << "\n~~~ Placement (best hit) Stats ~~~\n";
  std::cout << "Pendant lengths: \n";
  std::cout << "\tmin:\t" << minmax.min << "\n";
  std::cout << "\tmax:\t" << minmax.max << "\n";
  std::cout << "\tmedian:\t" << med << "\n";
  std::cout << "\tmean:\t" << meanstddev.mean << "\n";
  std::cout << "\tstddev:\t" << meanstddev.stddev << "\n";

  size_t discared_stddev_2 = 0;
  size_t discared_stddev_3 = 0;
  size_t discared_stddev_4 = 0;
  size_t discared_stddev_5 = 0;

  for( auto& sample : sampleset ) {
    for( auto& pq : sample ) {
      auto& p             = pq.placement_at( 0 );
      auto const edge_num = p.edge_num();
      auto z_score        = ( p.pendant_length - meanstddev.mean ) / meanstddev.stddev;

      if( z_score > 2 ) {
        discared_stddev_2++;
      }
      if( z_score > 3 ) {
        discared_stddev_3++;
      }
      if( z_score > 4 ) {
        discared_stddev_4++;
      }
      if( z_score > 5 ) {
        discared_stddev_5++;
        blame[ edge_num ].for_stddev_5++;
      }
    }
  }

  std::cout << "\n~~~ Outlier/Weirdo Detection ~~~\n";
  std::cout << "Best hits with pendant length more than " << thresh_mult << "x greater than:\n";
  std::cout << "\tLocal max:\t\t" << std::setprecision( 3 )
            << ( discarded_overmax / static_cast< double >( total_entries ) ) * 100 << "%\n";
  std::cout << "\tLocal average:\t\t" << std::setprecision( 3 )
            << ( discarded_overavg / static_cast< double >( total_entries ) ) * 100 << "%\n";
  std::cout << "\tLocal weighted average:\t" << std::setprecision( 3 )
            << ( discarded_overavg_weighted / static_cast< double >( total_entries ) ) * 100 << "%\n";

  std::cout << "Best hits with pendant length more than \n";
  std::cout << "\t2"
            << " sigma from the mean:\t" << std::setprecision( 3 )
            << ( discared_stddev_2 / static_cast< double >( total_entries ) ) * 100 << "%\n";
  std::cout << "\t3"
            << " sigma from the mean:\t" << std::setprecision( 3 )
            << ( discared_stddev_3 / static_cast< double >( total_entries ) ) * 100 << "%\n";
  std::cout << "\t4"
            << " sigma from the mean:\t" << std::setprecision( 3 )
            << ( discared_stddev_4 / static_cast< double >( total_entries ) ) * 100 << "%\n";
  std::cout << "\t5"
            << " sigma from the mean:\t" << std::setprecision( 3 )
            << ( discared_stddev_5 / static_cast< double >( total_entries ) ) * 100 << "%\n";

  std::cout << "\nEdges to blame for having placements with pendant length\n";
  std::cout << "\tover the local maximum:\n";
  std::sort( std::begin( blame ), std::end( blame ),
             []( blamestruct const& lhs, blamestruct const& rhs ) {
               return lhs.for_overmax > rhs.for_overmax;
             } );
  for( size_t i = 0; i < 5; ++i ) {
    std::cout << "\t\t" << std::setprecision( 2 )
              << ( blame[ i ].for_overmax / static_cast< double >( discarded_overmax ) ) * 100
              << "%\t" << blame[ i ].edge_num;
    if( not blame[ i ].name.empty() ) {
      std::cout << "\t" << blame[ i ].name;
    }
    std::cout << "\n";
  }

  std::cout << "\tover the local average:\n";
  std::sort( std::begin( blame ), std::end( blame ),
             []( blamestruct const& lhs, blamestruct const& rhs ) {
               return lhs.for_overavg > rhs.for_overavg;
             } );
  for( size_t i = 0; i < 5; ++i ) {
    std::cout << "\t\t" << std::setprecision( 2 )
              << ( blame[ i ].for_overavg / static_cast< double >( discarded_overavg ) ) * 100
              << "%\t" << blame[ i ].edge_num;
    if( not blame[ i ].name.empty() ) {
      std::cout << "\t" << blame[ i ].name;
    }
    std::cout << "\n";
  }

  std::cout << "\tover the local weighted average:\n";
  std::sort( std::begin( blame ), std::end( blame ),
             []( blamestruct const& lhs, blamestruct const& rhs ) {
               return lhs.for_overavg_weighted > rhs.for_overavg_weighted;
             } );
  for( size_t i = 0; i < 5; ++i ) {
    std::cout << "\t\t" << std::setprecision( 2 )
              << ( blame[ i ].for_overavg_weighted / static_cast< double >( discarded_overavg_weighted ) ) * 100
              << "%\t" << blame[ i ].edge_num;
    if( not blame[ i ].name.empty() ) {
      std::cout << "\t" << blame[ i ].name;
    }
    std::cout << "\n";
  }

  std::cout << "\tover 5 sigma from the mean:\n";
  std::sort( std::begin( blame ), std::end( blame ),
             []( blamestruct const& lhs, blamestruct const& rhs ) {
               return lhs.for_stddev_5 > rhs.for_stddev_5;
             } );
  for( size_t i = 0; i < 5; ++i ) {
    std::cout << "\t\t" << std::setprecision( 2 )
              << ( blame[ i ].for_stddev_5 / static_cast< double >( discared_stddev_5 ) ) * 100
              << "%\t" << blame[ i ].edge_num;
    if( not blame[ i ].name.empty() ) {
      std::cout << "\t" << blame[ i ].name;
    }
    std::cout << "\n";
  }

  // MatrixWriter<double>().to_stream( mat, std::cout, {}, {"lwr", "pendant_length"} );

  return 0;
}

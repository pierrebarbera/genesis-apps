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
#include <fstream>
#include <string>

#ifdef GENESIS_OPENMP
#include <omp.h>
#endif

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

void normalize( Sample& sample )
{
  auto const total = total_multiplicity( sample );

  std::for_each( std::begin( sample ),
                 std::end( sample ),
                 [total]( Pquery& pq ) {
                   for( auto& name : pq.names() ) {
                     name.multiplicity /= total;
                   }
                 } );
}

double medianPerqueryKRD( Sample const& lhs, Sample const& rhs )
{
  // normalize LWRs (to account for different LWR calc in raxml)
  // normalize_weight_ratios(lhs);
  // normalize_weight_ratios(rhs);

  // Collect the EMD distances and other parameters to make statistics about the results.
  std::vector< double > emd_results;
  std::vector< std::string > names_of_invalid_pqueries;

  // Prepare samples for the emd calculation by specifiying the underlying trees
  auto emd_lhs = Sample( lhs.tree() );
  auto emd_rhs = Sample( rhs.tree() );

  // For speedup, create maps from names of the samples to its pquery indices.
  auto name_map_l = std::unordered_map< std::string, size_t >();
  auto name_map_r = std::unordered_map< std::string, size_t >();

  for( size_t i = 0; i < lhs.size(); ++i ) {
    for( auto const& name_l : lhs.at( i ).names() ) {
      name_map_l[ name_l.name ] = i;
    }
  }

  for( size_t i = 0; i < rhs.size(); ++i ) {
    for( auto const& name_r : rhs.at( i ).names() ) {
      name_map_r[ name_r.name ] = i;
    }
  }

  // Iterate all pqueries of the left sample and find the equivalent pqueries of the right sample.
  for( auto& pqry_l : lhs ) {
    for( auto const& name_l : pqry_l.names() ) {

      // Check whether the right sample has a pquery with that name, and get it.
      if( name_map_r.count( name_l.name ) == 0 ) {
        continue;
      }
      auto& pqry_r = rhs.at( name_map_r[ name_l.name ] );

      // Add the pqueries to the emd samples.
      emd_lhs.add( pqry_l );
      emd_rhs.add( pqry_r );

      // Calculate the emd.
      emd_results.push_back( earth_movers_distance( emd_lhs, emd_rhs ) );

      // Remove the pqueries again.
      emd_lhs.clear_pqueries();
      emd_rhs.clear_pqueries();
    }
  }

  size_t n = emd_results.size() / 2;
  std::nth_element( emd_results.begin(), emd_results.begin() + n, emd_results.end() );
  return emd_results[ n ];
}

/**
 *  Outputs a pairwise median Phylogenetic Pantorovic-Rubinstein distance matrix for an arbitrary number of jplace files
 */
int main( int argc, char** argv )
{
  if( argc < 3 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + " <jplace-files...>" );
  }

  // In out dirs.
  std::vector< std::string > jplace_paths;
  for( int i = 1; i < argc; ++i ) {
    jplace_paths.push_back( std::string( argv[ i ] ) );
  }

  JplaceReader jplace_reader;
  auto sample_set = jplace_reader.read( from_files( jplace_paths ) );

  // get the average of the two trees to ensure comparability (this also readjusts the placement lengths)
  adjust_to_average_branch_lengths( sample_set );

  std::vector< std::string > names;
  for( auto const& s : sample_set.names() ) {
    names.push_back( s );
  }

  size_t n = sample_set.size();
  std::vector< std::pair< size_t, size_t > > idx;
  for( size_t i = 0; i < ( n - 1 ); ++i ) {
    for( size_t j = i + 1; j < n; ++j ) {
      idx.emplace_back( i, j );
    }
  }

  Matrix< double > krd_matrix( n, n, 0.0 );
#pragma omp parallel for schedule( dynamic )
  for( size_t k = 0; k < idx.size(); ++k ) {
    size_t i              = idx[ k ].first;
    size_t j              = idx[ k ].second;
    double krd            = medianPerqueryKRD( sample_set.at( i ), sample_set.at( j ) );
    krd_matrix.at( i, j ) = krd_matrix.at( j, i ) = krd;
  }

  MatrixWriter< double >().to_stream( krd_matrix, std::cout, {}, names );

  return 0;
}

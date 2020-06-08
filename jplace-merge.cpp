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

Sample read_and_merge( std::vector< std::string > const& paths )
{
  Sample result;
  size_t fc = 0;

// Read all jplace files and accumulate their pqueries.
#pragma omp parallel for schedule( dynamic )
  for( size_t fi = 0; fi < paths.size(); ++fi ) {
    auto const& cur = paths.at( fi );
    // User output.
    if( true ) {
#pragma omp critical( GAPPA_JPLACE_INPUT_PROGRESS )
      {
        ++fc;
        LOG_INFO << "Reading file " << fc << " of " << paths.size();
        LOG_INFO << ": " << cur << "\n";
      }
    }

    // Read in file. This is the part that can trivially be done in parallel.
    auto smpl = JplaceReader().read( from_file( cur ) );

    if( true ) {
      // normalize per-sample
      normalize( smpl );
    }

// The main merging is single threaded.
#pragma omp critical( GAPPA_JPLACE_INPUT_ACCUMULATE )
    {
      // Merge
      if( result.empty() ) {
        result = std::move( smpl );
      } else {
        try {
          // The function only throws if something is wrong with the trees.
          copy_pqueries( smpl, result );
        } catch( ... ) {
          throw std::runtime_error( "Input jplace files have differing reference trees." );
        }
      }
    }
  }

  merge_duplicates( result );

  if( true ) {
    normalize( result );
  }

  return result;
}

/**
 *  splits a jplace file into its  constituent samples, based on a standard OTU table
 */
int main( int argc, char** argv )
{

  utils::Options::get().number_of_threads( 4 );

  // Check if the command line contains the right number of arguments.
  if( argc < 2 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + " <jplace-files...>\n" );
  }

  // In out dirs.
  std::vector< std::string > jplace_paths;
  for( int i = 1; i < argc; ++i ) {
    jplace_paths.push_back( std::string( argv[ i ] ) );
  }

  JplaceReader jplace_reader;
  auto sample = read_and_merge( jplace_paths );

  JplaceWriter().write( sample, to_stream( std::cout ) );

  return 0;
}

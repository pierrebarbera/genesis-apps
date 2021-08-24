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

// #include <algorithm>
// #include <functional>
// #include <iomanip>
// #include <limits>
// #include <tuple>
// #include <utility>
#include <cmath>

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

/**
 * Remove PQueries that have any nan values 
 */
int main( int argc, char** argv )
{
  // Check if the command line contains the right number of arguments.
  if( argc < 2 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + " <jplace-files...> \n" );
  }

  std::vector< std::string > jplace_files;
  for( int i = 1; i < argc; ++i ) {
    jplace_files.emplace_back( argv[ i ] );
  }

  size_t total_trash = 0;

  for( auto& filename : jplace_files ) {
    size_t trash_count = 0;

    std::cout << "Cleaning " << filename << std::endl;
    Sample sample = JplaceReader().read( from_file( filename ) );
    Sample filtered_sample( sample.tree() );

    for( auto& pq : sample ) {

      Pquery filtered_pquery;
      for( auto& p : pq.placements() ) {

        if( std::isnan(p.like_weight_ratio) or
            std::isnan(p.likelihood) or
            std::isnan(p.proximal_length) or
            std::isnan(p.pendant_length)
         ) {
          ++trash_count;
        } else {
          filtered_pquery.add_placement( p );
        }
      }
      if( filtered_pquery.placement_size() ) {
        filtered_sample.add( filtered_pquery );
      }
    }

    auto filtered_name = "cleaned_" + filename;

    JplaceWriter().write( filtered_sample, to_file( filtered_name ) );

    std::cout << "Done! " << trash_count << " placements removed. Output: " << filtered_name << std::endl;

    total_trash += trash_count;
  }

  std::cout << "All done! Total removed placements: " << total_trash << std::endl;

  return 0;
}

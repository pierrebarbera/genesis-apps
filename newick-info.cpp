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
#include <vector>
#include <algorithm>

using namespace genesis;
using namespace genesis::sequence;
using namespace genesis::tree;
using namespace genesis::utils;

void fail( bool const fail_if_true, std::string const& msg )
{
  if( fail_if_true ) {
    throw std::runtime_error{ msg };
  }
}

int main( int argc, char** argv )
{
  // Check if the command line contains the right number of arguments.
  if( argc < 2 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + " <newick>..." );
  }

  Options::get().allow_file_overwriting( true );
  Logging::log_to_stdout();

  size_t const num_trees = argc - 1;

  double allsum_blavg = 0.0;
  double allsum_blsum = 0.0;
  double allsum_diam  = 0.0;
  double allsum_heigh = 0.0;

  for( int i = 1; i < argc; ++i ) {
    std::string tree_file( argv[ i ] );

    auto tree = CommonTreeNewickReader().read( from_file( tree_file ) );

    std::cout << "File: " << tree_file << "\n";
    std::cout << "  Tree is " << (is_bifurcating( tree ) ? "bifurcating" 
                                                        : "multifurcating") << "\n";
    std::cout << "  Tree is " << (is_rooted( tree ) ? "rooted" 
                                                   : "unrooted") << "\n";
    std::cout << "  Topology numbers:" << "\n";
    std::cout << "    leafs: " << std::to_string( leaf_node_count( tree ) ) << "\n";
    std::cout << "    inner: " << std::to_string( inner_node_count( tree ) ) << "\n";
    std::cout << "    edges: " << std::to_string( edge_count( tree ) ) << "\n";

    std::cout << "  Branch length stats, tree" << "\n";
    double blsum = length( tree );
    std::cout << "    BL sum:   " << std::to_string( blsum ) << "\n";
    double blavg = blsum / edge_count( tree );
    std::cout << "    BL avg:   " << std::to_string( blavg ) << "\n";
    double heigh = height( tree );
    std::cout << "    height:   " << std::to_string( heigh ) << "\n";
    double diam = diameter( tree );
    std::cout << "    diameter: " << std::to_string( diam ) << "\n";

    allsum_blavg += blavg;
    allsum_blsum += blsum;
    allsum_diam  += diam;
    allsum_heigh += heigh;
  }

  std::cout << "\n";

  std::cout << "Averages over all trees:" << "\n";
  std::cout << "  BL sum:   " << std::to_string( allsum_blsum / num_trees ) << "\n";
  std::cout << "  BL avg:   " << std::to_string( allsum_blavg / num_trees ) << "\n";
  std::cout << "  height:   " << std::to_string( allsum_heigh / num_trees ) << "\n";
  std::cout << "  diameter: " << std::to_string( allsum_diam / num_trees ) << "\n";
  return 0;
}

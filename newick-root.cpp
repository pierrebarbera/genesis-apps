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
#include <regex>

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

std::vector<TreeNode const*> nodes_by_pattern( Tree const& tree, std::string const& pattern )
{
  auto names = node_names( tree, true );
  std::regex reg( pattern );

  std::vector< std::string > root_names;

  for( auto const& name : names ) {
    std::smatch match;
    if( std::regex_search(  name, match, reg ) ) {
      std::cout << match.str() << std::endl;
      root_names.push_back( match.str() );
    }
  }

  return find_nodes( tree, root_names, true);
}

/**
 * LCA of given nodes
 */
TreeEdge& lowest_common_ancestor( Tree& tree, std::vector<TreeNode const*>& nodes )
{
    assert( not nodes.empty() );

    auto bipart = find_smallest_subtree( tree, bipartition_set( tree ), nodes );

    if ( bipart.empty() ) {
        throw std::invalid_argument{"Rooting could not be determined."};
    }

    return const_cast<TreeEdge&>( bipart.link().edge() );

}

int main( int argc, char** argv )
{
  // Check if the command line contains the right number of arguments.
  if( argc != 3 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + " <newick> <label-regex>" );
  }

  std::string tree_file( argv[ 1 ] );
  std::string root_pattern( argv[ 2 ] );

  auto tree = CommonTreeNewickReader().read( from_file( tree_file ) );

  // search for the leafs
  auto nodes = nodes_by_pattern( tree, root_pattern );
  // get lowest common ancestor edge of all found leafs
  auto& lca_branch = lowest_common_ancestor( tree, nodes );
  // root the tree there
  make_rooted( tree, lca_branch );

  auto writer = CommonTreeNewickWriter();
  writer.branch_length_precision(15);

  writer.write( tree, to_stream( std::cout ) );
  // std::cout << std::endl;

  return 0;
}

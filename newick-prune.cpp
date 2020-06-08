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

void fail( bool const fail_if_true, std::string const& msg )
{
  if( fail_if_true ) {
    throw std::runtime_error{ msg };
  }
}

int main( int argc, char** argv )
{
  // Check if the command line contains the right number of arguments.
  if( argc != 3 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + " <newick> <label>" );
  }

  Options::get().allow_file_overwriting( true );
  Logging::log_to_stdout();

  std::string tree_file( argv[ 1 ] );
  std::string prune_label( argv[ 2 ] );

  auto tree = CommonTreeNewickReader().read( from_file( tree_file ) );

  // search for the leaf node of the given label
  auto leaf_node =  find_node(tree, prune_label);
  auto& base_node = leaf_node->link().outer().node();

  fail( not leaf_node, "Given label was not found in the tree." );

  // remove leaf from tree
  delete_leaf_node( tree, *leaf_node );

  // remove inner node
  delete_linear_node( tree, base_node, []( TreeEdge& r_edge, TreeEdge& d_edge ) {
    r_edge.data< CommonEdgeData >().branch_length += d_edge.data< CommonEdgeData >().branch_length;
  } );

  auto writer = CommonTreeNewickWriter();
  writer.branch_length_precision(15);

  writer.to_stream( tree, std::cout );
  std::cout << std::endl;

  return 0;
}

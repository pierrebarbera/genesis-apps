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
#include <random>
#include <chrono>
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

/**
 * Checks if two trees have the same topology regardless of the order of their
 * links.
 * 
 * @param  lhs [description]
 * @param  rhs [description]
 * @return     [description]
 */
void validate( Tree const& lhs, Tree const& rhs )
{

  fail( lhs.node_count() != rhs.node_count(), "inequal node_count" );
  fail( lhs.link_count() != rhs.link_count(), "inequal link_count" );
  fail( lhs.edge_count() != rhs.edge_count(), "inequal edge_count" );

  // get names of all leafs
  auto lhs_names = node_names(lhs, true);
  auto rhs_names = node_names(rhs, true);

  fail( lhs_names.size() != rhs_names.size(), "inequal number of names" );

  // sort them for comparability
  std::sort( lhs_names.begin(), lhs_names.end() );
  std::sort( rhs_names.begin(), rhs_names.end() );

  // check that the names are equal
  for( size_t i = 0; i < lhs_names.size(); ++i ) {
    fail( lhs_names[i] != rhs_names[i], "unequal names" );
  }

  auto lhs_leafs =  find_nodes( lhs, lhs_names, true );
  auto rhs_leafs =  find_nodes( rhs, lhs_names, true );

  fail( lhs_leafs.size() != rhs_leafs.size(), "inequal node number" );

  // for each pair of leafs
  for( size_t i = 0; i < lhs_leafs.size(); ++i ) {
    for( size_t j = i; j < lhs_leafs.size(); ++j ) {
      if( i == j ) {
        continue;
      }

      // validate that the paths of both trees are the same

      auto it_l = path(*lhs_leafs[i], *lhs_leafs[j]).begin();
      auto it_r = path(*rhs_leafs[i], *rhs_leafs[j]).begin();

      auto it_l_end = path(*lhs_leafs[i], *lhs_leafs[j]).end();
      auto it_r_end = path(*rhs_leafs[i], *rhs_leafs[j]).end();

      for( ; it_l != it_l_end && it_r != it_l_end; ++it_l, ++it_r
      ) {
        fail( it_l.node().index() != it_r.node().index(),
              "indexcheck fail" );   
      }
    }
  }

}


int main( int argc, char** argv )
{
  // Check if the command line contains the right number of arguments.
  if( argc != 3 and argc != 2 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + " <newick> <optional seed>" );
  }

  Options::get().allow_file_overwriting( true );
  Logging::log_to_stdout();

  std::string tree_file( argv[ 1 ] );
  unsigned seed;
  if( argc == 3 ) {
    seed = std::stoi( argv[ 2 ] );
  } else {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
  }


  auto tree = CommonTreeNewickReader().read( from_file( tree_file ) );

  auto old_tree = tree;

  fail( not is_bifurcating( tree ), "Input tree must be strictly bifurcating!" );

  // init the random number generator
  std::default_random_engine gen( seed );
  std::uniform_int_distribution<int> dist(0,1);
  auto flip = std::bind(dist, gen);

  // make a list of internal nodes
  std::vector<TreeLink*> links_to_flip;
  for( auto it : postorder( tree ) ) {
    TreeLink* nlink   = &it.node().link();
    if( is_leaf( *nlink ) or it.is_last_iteration() ) {
      continue;
    }
    // select internal links whos children are to be swapped
    if( flip() ) {
      links_to_flip.push_back( nlink );
    }
  }

  // flip subtrees 
  for( auto nlink : links_to_flip ) {
    TreeLink* first   = &nlink->next().next();
    TreeLink* second  = &nlink->next();
    nlink->reset_next( first );
    first->reset_next( second );
    second->reset_next( nlink );
  }

  fail( not validate_topology( tree ), "Topology invalid." );

  #ifdef DDEBUG
  validate( old_tree, tree );
  #endif
  
  auto writer = CommonTreeNewickWriter();
  writer.branch_length_precision(20);
  writer.write( tree, to_stream( std::cout ) );
  std::cout << std::endl;

  return 0;
}

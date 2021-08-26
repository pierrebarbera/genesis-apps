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

#include "common.hpp"

#include <algorithm>
#include <fstream>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

using namespace genesis;
using namespace genesis::sequence;
using namespace genesis::utils;

bool has( const std::vector< size_t >& cont, size_t val )
{
  return ( std::find( std::begin( cont ), std::end( cont ), val ) != std::end( cont ) );
}

std::vector< size_t > get_rand_unique( const size_t num, const size_t min, const size_t max )
{
  std::vector< size_t > ret;
  ret.reserve( num );

  assert( min <= max );

  std::random_device rd;
  std::mt19937 dom( rd() );
  std::uniform_int_distribution< size_t > ran( min, max );

  size_t to_add = num;
  while( to_add > 0 ) {
    auto rand_num = ran( dom );
    if( not has( ret, rand_num ) ) {
      ret.push_back( rand_num );
      --to_add;
    }
  }

  return ret;
}

void remove_duplicates( SequenceSet& set )
{
  std::unordered_set< std::string > labels;

  auto new_end = std::remove_if(
      std::begin( set ),
      std::end( set ),
      [ & ]( Sequence const& seq ) {
        return ( not labels.insert( seq.label() ).second );
      } );

  set.remove( new_end, std::end( set ) );
}

int main( int argc, char** argv )
{
  // Check if the command line contains the right number of arguments.
  if( argc < 2 or argc > 3 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + " <number of sequences> <seqfile|stdin>" );
  }

  // Get labels of reference alignment.
  auto in_set = read_any_seqfile( argv[ 2 ], (argc == 2) );

  remove_duplicates( in_set );

  const auto num = static_cast< size_t >( std::stoi( argv[ 1 ] ) );

  // get random indices
  auto idx = get_rand_unique( num, 0, in_set.size() );

  FastaOutputIterator out { to_stream( std::cout ) };
  for( auto i : idx ) {
    out << in_set[ i ];
  }

  return 0;
}

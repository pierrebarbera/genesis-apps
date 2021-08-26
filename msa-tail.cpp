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

#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace genesis;
using namespace genesis::sequence;
using namespace genesis::utils;

int main( int argc, char** argv )
{
  // Check if the command line contains the right number of arguments.
  if( argc < 2 or argc > 3 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + " <n> <fasta_msa|stdin>" );
  }

  SequenceSet in_set;
  FastaReader().read( (argc == 2) ? from_stream( std::cin ) : from_file( argv[ 2 ] ), in_set );

  bool const first_n = (argv[ 1 ][0] == '+');

  auto const n_in = std::stoi( argv[ 1 ] );

  if( n_in <= 0 ) {
    throw std::runtime_error( "n must be nonnegative" );
  }

  auto const n = static_cast<size_t>(n_in);

  auto const max  = std::min( n, in_set.size() );
  auto const skip = (first_n) ? n : in_set.size() - max;

  FastaOutputIterator out { to_stream( std::cout ) };
  for( size_t i = skip; i < in_set.size(); ++i ) {
    out << in_set[ i ];
  }

  return 0;
}

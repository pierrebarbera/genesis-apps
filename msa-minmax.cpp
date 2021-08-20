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

#include <iostream>
#include <string>

using namespace genesis;
using namespace genesis::sequence;
using namespace genesis::utils;

int main( int argc, char** argv )
{
  // Check if the command line contains the right number of arguments.
  if( argc != 4 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + " <fasta_msa> <min-N> <max-N>" );
  }

  auto const infile = std::string( argv[ 1 ] );
  auto const min_N  = std::stoi( argv[ 2 ] );
  auto const max_N  = std::stoi( argv[ 3 ] );

  if( min_N < 0 ) {
    throw std::runtime_error("min-N must be positive");
  }
  if( max_N < 0 ) {
    throw std::runtime_error("max-N must be positive");
  }

  auto fasta_in = FastaInputIterator( from_file( infile ) );
  FastaOutputIterator fasta_out { to_stream( std::cout ) };

  while( fasta_in ) {
    auto seq = *fasta_in;

    if( seq.sites().empty() || seq.label().empty() ) {
      throw std::runtime_error( "Invalid sequences with empty label or sites." );
    }

    // only return sequences that have between min and max nucleotides
    auto len = find_sites( seq, nucleic_acid_codes_all_letters() ).count();
    if( len >= min_N && len <= max_N ) {
      fasta_out << seq;
    }

    ++fasta_in;
  }

  return 0;
}

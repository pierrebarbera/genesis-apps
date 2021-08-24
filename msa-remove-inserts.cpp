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
#include <limits>

using namespace genesis;
using namespace genesis::sequence;
using namespace genesis::utils;

/*
  Drop insertion columns from alignment.
 */
int main( int argc, char** argv )
{
  // Check if the command line contains the right number of arguments.
  if( argc < 2 or argc > 3 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + " <fasta_msa> [<num-seqs>]" );
  }

  auto const infile     = std::string( argv[ 1 ] );
  auto const trimchars  = ".";
  size_t const num_seqs = (argc == 3) ? std::stoi( argv[ 2 ] ) : std::numeric_limits< size_t >::max();

  if (num_seqs == 0) {
    throw std::runtime_error("num_seqs of 0 makes no sense");
  }

  auto fasta_in = FastaInputIterator(
    (infile == "--") ? from_stream( std::cin ) : from_file( infile ) );

  if (not fasta_in) {
    throw std::runtime_error("input fasta empty/invalid?");
  }

  auto const msa_width = fasta_in->length();

  Bitvector insert_positions( msa_width, false );

  // go through the file and identify columns that correspond to inserts ("." if coming from hmmer)
  size_t cur_seq = 0;
  while( fasta_in and cur_seq < num_seqs ) {
    auto seq = *fasta_in;

    if( seq.sites().empty() || seq.label().empty() ) {
      throw std::runtime_error( "Invalid sequences with empty label or sites." );
    }

    if( seq.length() != msa_width ) {
      throw std::runtime_error( "Sequence had incorrect length (file not alignment?)" );
    }

    insert_positions |= find_sites( seq, trimchars );

    ++fasta_in;
    ++cur_seq;
  }

  // go through again, write only non-insert columns
  fasta_in = FastaInputIterator( from_file( infile ) );
  FastaOutputIterator fasta_out { to_stream( std::cout ) };
  while( fasta_in ) {
    auto seq = *fasta_in;

    remove_sites( seq, insert_positions );

    fasta_out << seq;
    
    ++fasta_in;
  }
  return 0;
}

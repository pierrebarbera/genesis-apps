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
  Drop columns columns from alignment.
 */
int main( int argc, char** argv )
{
  // Check if the command line contains the right number of arguments.
  if( argc < 2 or argc > 3 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + " <fasta_msa> [<first-n>]" );
  }

  auto const infile     = std::string( argv[ 1 ] );
  size_t const first_n  = (argc == 3) ? std::stoi( argv[ 2 ] ) : std::numeric_limits< size_t >::max();

  if (first_n == 0) {
    throw std::runtime_error("first_n of 0 makes no sense");
  }

  auto fasta_in = FastaInputIterator(
    (infile == "-") ? from_stream( std::cin ) : from_file( infile ) );

  if (not fasta_in) {
    throw std::runtime_error("input fasta empty/invalid?");
  }

  auto const msa_width = fasta_in->length();

  auto site_count_mat = SiteCounts( nucleic_acid_codes_all(), msa_width );

  // go through the (first n) sequences and identify char counts per site 
  size_t cur_seq = 0;
  while( fasta_in and cur_seq < first_n ) {
    auto seq = *fasta_in;

    if( seq.sites().empty() || seq.label().empty() ) {
      throw std::runtime_error( "Invalid sequences with empty label or sites." );
    }

    if( seq.length() != msa_width ) {
      throw std::runtime_error( "Sequence had incorrect length (file not alignment?)" );
    }

    site_count_mat.add_sequence( seq, false );

    ++fasta_in;
    ++cur_seq;
  }

  assert(site_count_mat.added_sequences_count() == first_n);

  // create a bitvector indicating which columns to remove (those that were static) 
  Bitvector to_remove( msa_width, false );
  for (size_t site = 0; site < msa_width; ++site) {
    // first determine if the current site is static:
    // this is the case one nt count equals the total number of sequences
    bool is_static = false;
    for( char const nt : site_count_mat.characters() ) {
      if( site_count_mat.count_of( nt, site ) == first_n ) {
        is_static = true;
        break;
      }
    }

    to_remove.set( site, is_static );
  }

  // go through again, write only non-static columns
  fasta_in = FastaInputIterator( from_file( infile ) );
  FastaOutputIterator fasta_out { to_stream( std::cout ) };
  while( fasta_in ) {
    auto seq = *fasta_in;

    remove_sites( seq, to_remove );

    fasta_out << seq;
    
    ++fasta_in;
  }
  return 0;
}

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
        std::string( "Usage: " ) + argv[ 0 ] + " <fasta_msa> <trim_left> <trim_right>" );
  }

  auto const infile     = std::string( argv[ 1 ] );
  auto const trim_left  = std::stoi( argv[ 2 ] );
  auto const trim_right = std::stoi( argv[ 3 ] );

  if( ( trim_left < 0 ) or ( trim_right < 0 ) ) {
    throw std::runtime_error( "Both trim values must be >= 0!" );
  }

  auto const total_trimmed = trim_left + trim_right;

  auto fasta_in = FastaInputIterator( from_file( infile ) );
  FastaOutputIterator fasta_out { std::cout };

  while( fasta_in ) {
    auto seq = *fasta_in;

    if( seq.sites().empty() || seq.label().empty() ) {
      throw std::runtime_error( "Invalid sequences with empty label or sites." );
    }

    if( total_trimmed >= seq.length() ) {
      throw std::runtime_error( "Trying to trim a sequence out of existence." );
    }

    auto begin = seq.sites().begin() + trim_left;
    auto end   = seq.sites().end() - trim_right;

    seq.sites( std::string( begin, end ) );

    fasta_out = seq;

    ++fasta_in;
  }

  return 0;
}

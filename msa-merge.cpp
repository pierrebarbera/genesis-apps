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

#include <fstream>
#include <string>

using namespace genesis;
using namespace genesis::sequence;
using namespace genesis::utils;

/**
 * Takes a set of fasta files and merges then, prepending the first part of the
 * filename to each sequence contained
 * 
 * @param  argc [description]
 * @param  argv [description]
 * @return      [description]
 */
int main( int argc, char** argv )
{
  // Check if the command line contains the right number of arguments.
  if( argc < 2 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + "msa_files..." );
  }

  FastaOutputIterator fasta_out { to_stream( std::cout ) };

  // read the files
  for( int i = 1; i < argc; ++i ) {
    std::string const file_path = argv[ i ];
    auto fasta_in = FastaInputIterator( from_file( file_path ) );

    auto filename = split( file_filename( file_basename( file_path ) ), ".");
    if ( filename.empty() ){
      throw std::runtime_error(
        std::string( "File didn't contain suitable prefix in its name: " ) 
        + file_path ); 
    }

    auto prefix = filename[0];

    while( fasta_in ) {
      auto seq = *fasta_in;

      // prefix the sequence
      seq.label( prefix + "." + seq.label() );

      fasta_out << seq;

      ++fasta_in;
    }
  }
  return 0;
}

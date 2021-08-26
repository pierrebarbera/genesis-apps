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
using namespace genesis::utils;

int main( int argc, char** argv )
{
  // Check if the command line contains the right number of arguments.
  if( argc < 2 or argc > 3 ) {
    LOG_INFO << "Usage: " << argv[ 0 ] << " <n> <fasta|stdin>";
    return 1;
  }

  const size_t num = std::stoi( argv[ 1 ] );
  auto seq_file    = std::string( argv[ 2 ] );
  auto writer      = FastaWriter();

  auto iter = FastaInputIterator( (argc == 2) ? from_stream( std::cin ) : from_file( seq_file ) );
  for( size_t i = 0; iter and ( i < num ); ++i ) {
    writer.write_sequence( *iter, std::cout );
    iter.increment();
  }

  return 0;
}

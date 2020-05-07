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

#include "common.hpp"

#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace genesis;
using namespace genesis::sequence;
using namespace genesis::utils;

int main( int argc, char** argv )
{
  // Activate logging.
  utils::Logging::log_to_stdout();
  LOG_INFO << "Started " << utils::current_time();

  // Check if the command line contains the right number of arguments.
  if( argc < 2 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + "msa_files..." );
  }

  for( int i = 1; i < argc; ++i ) {
    std::string const msa_file = argv[ i ];
    // Get labels of reference alignment.
    auto set   = read_any_seqfile( msa_file );
    bool valid = true;
    LOG_INFO << "File: " << msa_file;

    auto const& first_seq = set[ 0 ];

    auto sites = first_seq.size();
    LOG_INFO << "Sites: " << sites;

    // bitvector indicating for each site wether there was more than one type of nucleotide
    Bitvector site_is_variable( sites );

    size_t n = 0;
    for( auto& s : set ) {
      n++;
      if( s.size() != sites ) {
        LOG_ERR << "Incorrect number of sites for sequence " << n;
        LOG_ERR << sites << " vs " << s.size();
        valid = false;
        continue;
      }

      // check for site variability
      for( size_t k = 0; k < sites; ++k ) {
        if( s[ k ] != first_seq[ k ] ) {
          site_is_variable.set( k );
        }
      }
    }

    LOG_INFO << "Variable Sites: " << site_is_variable.count();

    LOG_INFO << "Sequences: " << n;
    if( valid ) {
      LOG_INFO << "File OK!";
    } else {
      LOG_INFO << "File NOT OK!";
    }
  }
  LOG_INFO << "Finished " << utils::current_time();
  return 0;
}

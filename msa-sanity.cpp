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

constexpr char VAR_ONCE = '1';
constexpr char VAR_MANY = '2';

struct counts {
  size_t variable  = 0;
  size_t singleton = 0;
  size_t constant  = 0;
};

counts eval_sites( std::string const& sites )
{
  counts ret;

  for(auto const& s : sites){
    switch(s) {
      case VAR_ONCE:
        ret.singleton++;
        __attribute__ ((fallthrough));
      case VAR_MANY:
        ret.variable++;
        break;
    }
  }

  ret.constant = sites.size() - ret.variable;

  return ret;
}

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
    LOG_INFO << " ======================= ";
    LOG_INFO << "File: " << msa_file;

    normalize_nucleic_acid_codes( set[ 0 ] );
    auto const& first_seq = set[ 0 ];

    auto sites = first_seq.size();

    // make an editable sequence where we set the sites to special chars if it is variable

    std::string check_sites = first_seq.sites();

    size_t n = 0;
    for( auto& s : set ) {
      // make sure the seqs all look the same
      normalize_nucleic_acid_codes( s );
      n++;
      if( s.size() != sites ) {
        LOG_ERR << "Incorrect number of sites for sequence " << n;
        LOG_ERR << sites << " vs " << s.size();
        valid = false;
        continue;
      }

      // check for site variability
      for( size_t k = 0; valid and k < sites; ++k ) {
        auto& check_site = check_sites[ k ];
        if( s[ k ] != check_site ) {
          check_site = ( check_site == VAR_ONCE ) ? VAR_MANY : VAR_ONCE;
        }
      }
    }

    auto count = eval_sites( check_sites );

    LOG_INFO << n << " sequences";

    LOG_INFO << sites << " sites";
    LOG_INFO << "\t" << count.variable << "\t" << "variable";
    LOG_INFO << "\t" << count.singleton << "\t" << "singleton";
    LOG_INFO << "\t" << count.constant << "\t" << "constant";
    LOG_INFO << "\t" << gap_sites( set ).count() << "\t" << "gap";

    if( valid ) {
      LOG_INFO << "File OK!";
    } else {
      LOG_INFO << "File NOT OK!";
    }
  }
  LOG_INFO << "Finished " << utils::current_time();
  return 0;
}

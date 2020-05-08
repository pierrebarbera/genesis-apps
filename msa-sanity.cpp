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

struct counts_t {
  size_t variable    = 0;
  size_t pars_inform = 0;
  size_t singleton   = 0;
  size_t constant    = 0;
};

counts_t eval_sites( SiteCounts const& counts )
{
  counts_t ret;

  auto const sites     = counts.length();
  auto const num_chars = counts.characters().size();

  for( size_t site_idx = 0; site_idx < sites; ++site_idx ) {
    size_t num_nonzero_chars = 0;
    size_t num_two_or_more   = 0;
    for( size_t char_idx = 0; char_idx < num_chars; ++char_idx ) {
      auto const c = counts.count_at( char_idx, site_idx );
      num_nonzero_chars += static_cast< bool >( c );
      num_two_or_more += ( c > 1 );
    }

    if( num_nonzero_chars > 1 ) {
      ret.variable++;
      if( num_two_or_more > 1 ) {
        ret.pars_inform++;
      } else if( num_two_or_more == 1 ) {
        ret.singleton++;
      }
    } else {
      ret.constant++;
    }
  }

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

    SiteCounts site_counts( nucleic_acid_codes_all(), sites );

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
      site_counts.add_sequence( s, false );
    }

    auto count = eval_sites( site_counts );

    LOG_INFO << n << " sequences";

    LOG_INFO << sites << " sites";
    LOG_INFO << "\t" << count.variable << "\t"
             << "variable";
    LOG_INFO << "\t" << count.pars_inform << "\t"
             << "parsimony-informative";
    LOG_INFO << "\t" << count.singleton << "\t"
             << "singleton";
    LOG_INFO << "\t" << count.constant << "\t"
             << "constant";
    LOG_INFO << "\t" << gap_sites( set ).count() << "\t"
             << "gap";

    if( valid ) {
      LOG_INFO << "File OK!";
    } else {
      LOG_INFO << "File NOT OK!";
    }
  }
  LOG_INFO << "Finished " << utils::current_time();
  return 0;
}

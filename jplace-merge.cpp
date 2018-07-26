/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2017 Lucas Czech

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
    Lucas Czech <lucas.czech@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "genesis/genesis.hpp"

#include <fstream>
#include <string>
#include <algorithm>

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

void normalize(Sample& sample) {
    auto const total = total_multiplicity(sample);

    std::for_each(  std::begin(sample),
                    std::end(sample),
                    [total](Pquery& pq){
                        for (auto& name : pq.names()) {
                            name.multiplicity /= total;
                        }
                    });
}

Sample read_and_merge( std::vector<std::string> const& paths )
{
    Sample result;
    size_t fc = 0;

    // Read all jplace files and accumulate their pqueries.
    #pragma omp parallel for schedule(dynamic)
    for( size_t fi = 0; fi < paths.size(); ++fi ) {
        auto const& cur = paths.at( fi );
        // User output.
        if( true ) {
            #pragma omp critical(GAPPA_JPLACE_INPUT_PROGRESS)
            {
                ++fc;
                LOG_INFO << "Reading file " << fc << " of " << paths.size();
                LOG_INFO << ": " << cur << "\n";
            }
        }

        // Read in file. This is the part that can trivially be done in parallel.
        auto smpl = JplaceReader().from_file( cur );

        if (true) {
            // normalize per-sample
            normalize(smpl);
        }

        // The main merging is single threaded.
        #pragma omp critical(GAPPA_JPLACE_INPUT_ACCUMULATE)
        {
            // Merge
            if( result.empty() ) {
                result = std::move( smpl );
            } else {
                try{
                    // The function only throws if something is wrong with the trees.
                    copy_pqueries( smpl, result );
                } catch( ... ) {
                    throw std::runtime_error( "Input jplace files have differing reference trees." );
                }
            }

        }
    }

    merge_duplicates( result );

    if (true) {
        normalize(result);
    }

    return result;
}

/**
 *  splits a jplace file into its  constituent samples, based on a standard OTU table
 */
int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.date = true;
    utils::Logging::details.time = true;

    utils::Options::get().number_of_threads( 4 );
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    LOG_INFO << "Started";

    // Check if the command line contains the right number of arguments.
    if ( argc < 3 ) {
        throw std::runtime_error(
            "Usage: jplace-merge <out-name> <jplace-files...>\n"
        );
    }

    // In out dirs.
    std::vector<std::string> jplace_paths;
    for (int i = 2; i < argc; ++i) {
        jplace_paths.push_back( std::string( argv[i] ) );
    }
    auto outfile = std::string( argv[1] ) + ".jplace";

    JplaceReader jplace_reader;
    auto sample = read_and_merge( jplace_paths );

    JplaceWriter().to_file( sample, outfile );

    LOG_INFO << "Finished";
    return 0;
}

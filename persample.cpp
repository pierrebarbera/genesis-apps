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

/**
 *  splits a jplace file into its  constituent samples, based on a standard OTU table
 */
int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;

    LOG_INFO << "Started";

    // Check if the command line contains the right number of arguments.
    if (argc != 4) {
        throw std::runtime_error(
            "Usage: persample <jplace-file> <otu-table> <out-path>\n"
        );
    }

    // In out dirs.
    auto jplacefile = std::string( argv[1] );
    auto otufile = std::string( argv[2] );
    auto outdir = utils::dir_normalize_path( std::string( argv[3] ) );
    utils::dir_create(outdir);

    // -------------------------------------------------------------------------
    //     demangle jplace file based on OTU table
    // -------------------------------------------------------------------------
    JplaceReader jplace_reader;
    auto in_sample = jplace_reader.from_file( jplacefile );

    LOG_INFO << "Finished reading input jplace file: " << jplacefile;

    CsvReader csvreader;

    csvreader.separator_chars("\t");
    csvreader.comment_chars("#");

    auto table = csvreader.from_file( otufile );

    LOG_INFO << "Finished reading OTU table: " << otufile;


    auto& headers = table[0];
    LOG_INFO << "Splitting into " << headers.size() - 1 << " separate sample files.";

    auto writer = JplaceWriter();
    // iterate over the OTU table, column major
    for (size_t col = 1; col < table[0].size() - 1; ++col) {
        // create an output sample
        auto sample_id = headers[col];
        Sample out_sample( in_sample.tree() );

        LOG_DBG << "Sample: " << sample_id ;

        for (size_t row = 1; row < table.size(); ++row) {
            // get the OTU-id
            auto otu_id         = table[row][0];
            auto multiplicity   = std::stod( table[row][col] );

            // only write if the otu was in the sample to begin with
            if ( multiplicity > 0 ) {
                // look up the relevant PQuery
                auto it = std::find_if( std::begin(in_sample), std::end(in_sample),
                    [otu_id]( const Pquery& pq ) {
                        bool ret = false;
                        for ( auto& name : pq.names() ) {
                            if ( otu_id == name ) {
                                ret = true;
                            }
                        }
                        return ret;
                    }

                );

                if ( it != std::end(in_sample) ) {
                    // add pquery to output sample including multiplicity count
                    auto& pq = out_sample.add( *it );
                    // we need to manually set name and multiplicity such that other names/multiplicities
                    // don't get copied over from the input sample
                    pq.clear_names();
                    pq.add_name( otu_id, multiplicity );
                }
            }
        }

        writer.to_file( out_sample, outdir + sample_id + ".jplace" );
    }

    LOG_INFO << "Finished";
    return 0;
}

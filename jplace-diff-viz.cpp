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
    utils::Logging::details.date = true;
    utils::Logging::details.time = true;

    utils::Options::get().number_of_threads( 4 );
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    LOG_INFO << "Started";

    // Check if the command line contains the right number of arguments.
    if ( argc != 3 ) {
        throw std::runtime_error(
            "Usage: persample <jplace-file> <jplace-file>\n"
        );
    }

    // In out dirs.
    std::vector<std::string> jplace_paths;
    for (int i = 2; i < argc; ++i) {
        jplace_paths.push_back( std::string( argv[i] ) );
    }
    auto outfile = std::string( argv[1] );
    // auto outdir = utils::dir_normalize_path( std::string( "." ));
    // utils::dir_create(outdir);


    JplaceReader jplace_reader;
    auto const sample_set = jplace_reader.from_files( jplace_paths );

    auto const pwdmat = earth_movers_distance(sample_set);

    std::ofstream out(outfile);
    out << pwdmat;

    LOG_INFO << "Finished";
    return 0;
}

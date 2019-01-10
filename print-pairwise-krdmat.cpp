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

#include <fstream>
#include <string>
#include <algorithm>

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

/**
 *  Outputs a pairwise Phylogenetic Kantorovic-Rubinstein distance matrix for an arbitrary number of jplace files
 */
int main( int argc, char** argv )
{
    if ( argc < 2 ) {
        throw std::runtime_error(
            std::string("Usage: ") + argv[0] +  " <jplace-files...>"
        );
    }

    // In out dirs.
    std::vector<std::string> jplace_paths;
    for (int i = 1; i < argc; ++i) {
        jplace_paths.push_back( std::string( argv[i] ) );
    }

    JplaceReader jplace_reader;
    auto sample_set = jplace_reader.read( from_files( jplace_paths ) );

    for (auto& sample : sample_set) {
        normalize(sample);
    }

    auto const pwdmat = earth_movers_distance(sample_set);

    std::vector<std::string> names;
    for ( auto const& name : sample_set.names() ) {
        names.push_back( name );
    }

    MatrixWriter<double>().to_stream( pwdmat, std::cout, {}, names );

    return 0;
}

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
using namespace genesis::tree;
using namespace genesis::utils;

int main( int argc, char** argv )
{
    // Check if the command line contains the right number of arguments.
    if (argc != 3) {
        throw std::runtime_error(
            std::string( "Usage: " ) + argv[0] + " <type> <file>"
        );
    }

    std::string filetype(argv[1]);
    std::string filename(argv[2]);

    if (to_lower(filetype) == "fasta") {
        // Prepare reading and writing files.
        SequenceSet set;

        // Get labels of reference alignment.
        FastaReader().from_file(filename , set );

        for (auto& seq : set) {
            seq.label( SHA1().from_string_hex( seq.label() ) );
        }

        FastaWriter().to_stream(set, std::cout);
    } else if (to_lower(filetype) == "newick") {
        // Get labels of reference alignment.
        auto tree = DefaultTreeNewickReader().from_file( filename );

        auto leaf_ids = leaf_node_indices( tree );

        for (auto id : leaf_ids) {
            auto& name = tree.node_at( id ).data<DefaultNodeData>().name;
            name = SHA1().from_string_hex( name );
        }

        DefaultTreeNewickWriter().to_stream(tree, std::cout);
    }

    return 0;
}

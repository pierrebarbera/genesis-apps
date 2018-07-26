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
#include <unordered_map>
#include <vector>
#include <random>
#include <algorithm>

using namespace genesis;
using namespace genesis::sequence;

bool has(const std::vector<size_t>& cont, size_t val)
{
    return ( std::find(std::begin(cont), std::end(cont), val) != std::end(cont) );
}

std::vector<size_t> get_rand_unique(const size_t num, const size_t min, const size_t max)
{
    std::vector<size_t> ret;
    ret.reserve(num);

    assert(min <= max);

    std::random_device rd;
    std::mt19937 dom(rd());
    std::uniform_int_distribution<size_t> ran(min, max);

    size_t to_add = num;
    while (to_add > 0) {
        auto rand_num = ran(dom);
        if ( not has(ret, rand_num) ) {
            ret.push_back(rand_num);
            --to_add;
        }
    }

    return ret;
}

void remove_duplicates(SequenceSet& set)
{
    std::unordered_set<std::string> labels;

    auto new_end =
        std::remove_if(
            std::begin(set),
            std::end(set),
            [&] ( Sequence const& seq ) {
                return ( not labels.insert(seq.label()).second );
            }
        );

    set.remove(new_end, std::end(set));
}

int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    // utils::Logging::log_to_stdout();
    // LOG_INFO << "Started " << utils::current_time();

    // Check if the command line contains the right number of arguments.
    if (argc != 3) {
        throw std::runtime_error(
            std::string( "Usage: " ) + argv[0]  + " fasta_msa" + " <number of sequences>" 
        );
    }

    // Prepare reading and writing files.
    auto reader = FastaReader();
    auto writer = FastaWriter();
    auto in_set = SequenceSet();
    auto out_set = SequenceSet();

    // Get labels of reference alignment.
    reader.from_file( argv[1], in_set );

    remove_duplicates(in_set);

    const auto num = static_cast<size_t>(std::stoi(argv[2]));

    // get random indices
    auto idx = get_rand_unique(num, 0, in_set.size());

    for (auto i : idx) {
        out_set.add(in_set[i]);
    }

    writer.to_stream(out_set, std::cout);

    // LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
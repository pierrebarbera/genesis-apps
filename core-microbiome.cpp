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

struct CoreEntry {
    CoreEntry(double i, double fract) : incidence(i), significance(1, fract) {};
    double incidence;
    std::vector<double> significance;
};

using Core = std::unordered_map<std::string, CoreEntry>;

void parse(Core& lhs, CsvReader::Table& rhs) {
    bool first = true;
    for (auto& row : rhs) {
        // skip header
        if (first) {
            first=false;
            continue;
        }
        auto& taxopath = row.at(4);
        double fract = atof( row.at(1).c_str() );
        auto it = lhs.find( taxopath );

        if ( it == lhs.end() ) {
            lhs.emplace(taxopath, CoreEntry(1.0, fract) );
        } else {
            it->second.incidence += 1.0;
            it->second.significance.emplace_back( fract );
        }
    }
}

void to_percent_and_filter( Core& core, size_t const total, double const thresh ) {
    std::vector<std::string> erase_keys;
    for (auto i = core.begin(); i != core.end(); ++i) {
        i->second.incidence /= static_cast<double>(total);
        if (i->second.incidence < thresh) {
            erase_keys.emplace_back( i->first );
        }
    }

    // erasure phase
    for (auto& key : erase_keys) {
        core.erase( key );
    }
}

void to_stream(Core& core, std::ostream& stream) {
    stream << "taxopath\tincidence\tsignificance\n";
    for (auto i = core.begin(); i != core.end(); ++i) {
        stream << i->first << "\t";
        stream << i->second.incidence << "\t";

        auto& sigvec = i->second.significance;
        double total = std::accumulate(sigvec.begin(), sigvec.end(), 0.0);
        double avg_sig = total / (double)sigvec.size();
        stream << avg_sig << "\n";
    }
}

/**
 *  determines a core profile of two or more taxonomic profiles as generated with gappa analyze assign
 */
int main( int argc, char** argv )
{
    // Check if the command line contains the right number of arguments.
    if ( argc < 3 ) {
        throw std::runtime_error(
            "Usage: core-microbiome <inclusion-threshold> <jplace-files...>\n"
            "Where inclusion threshold is the percentage [1.0,0.0) of samples a taxon has to appear in to be included"
        );
    }

    // In out dirs.
    std::vector<std::string> profile_files;
    for (int i = 2; i < argc; ++i) {
        profile_files.push_back( std::string( argv[i] ) );
    }
    auto thresh = std::stod( argv[1] );

    // read in the profiles
    CsvReader reader;
    reader.separator_chars("\t");
    std::vector<CsvReader::Table> profiles;
    for (auto const& p : profile_files) {
        profiles.emplace_back( reader.from_file( p ) );
    }

    Core result;
    for (auto& profile : profiles) {
        parse(result, profile);
    }

    to_percent_and_filter( result, profiles.size(), thresh );

    to_stream( result, std::cout );

    return 0;
}

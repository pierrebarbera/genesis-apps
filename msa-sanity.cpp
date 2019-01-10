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
#include <unordered_map>
#include <vector>

using namespace genesis;
using namespace genesis::sequence;
using namespace genesis::utils;

int main( int argc, char** argv )
{
    bool valid = true;

    // Activate logging.
    utils::Logging::log_to_stdout();
    LOG_INFO << "Started " << utils::current_time();

    // Check if the command line contains the right number of arguments.
    if (argc != 2) {
        throw std::runtime_error(
            std::string( "Usage: " ) + argv[0] + "fasta_msa"
        );
    }


    // Prepare reading and writing files.
    auto reader = FastaReader();
    auto set = SequenceSet();

    // Get labels of reference alignment.
    reader.read( from_file( argv[1] ), set );

    auto sites = set[0].size();
    LOG_INFO << "Sites: " << sites;
    size_t i = 0;
    for (auto& s : set)
    {
        if (s.size() != sites)
        {
            LOG_ERR << "Incorrect number of sites for sequence " << i;
            LOG_ERR << sites << " vs " << s.size();
            valid=false;
        }
        i++;
    }
    LOG_INFO << "Sequences: " << i;
    if (valid) LOG_INFO << "File OK!"; else LOG_INFO << "File NOT OK!";
    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}

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

#include <algorithm>
#include <fstream>
#include <string>

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

void print_tree( Sample const& sample, std::string const& filename )
{
  auto layout         = LayoutParameters{};
  layout.stroke.width = 15.0;
  layout.ladderize    = false;

  auto edge_masses = placement_mass_per_edges_with_multiplicities( sample );
  auto cmap        = ColorMap( color_list_viridis() );
  auto cnorm       = ColorNormalizationLinear();
  cnorm.autoscale( edge_masses );
  auto const colors = cmap( cnorm, edge_masses );

  write_color_tree_to_svg_file( sample.tree(), layout, colors, cmap, cnorm, filename );
}

/**
 *  Makes a tree visualizing the differences between two jplace
 */
int main( int argc, char** argv )
{
  Options::get().allow_file_overwriting( true );
  // Activate logging.
  Logging::log_to_stdout();
  Logging::details.date = true;
  Logging::details.time = true;

  Options::get().number_of_threads( 4 );
  LOG_BOLD << utils::Options::get().info();
  LOG_BOLD;

  LOG_INFO << "Started";

  // Check if the command line contains the right number of arguments.
  if( argc != 4 ) {
    throw std::runtime_error(
        std::string( "Usage: " ) + argv[ 0 ] + " <jplace-file> <jplace-file> <out-name>\n" );
  }
  auto outfile = std::string( argv[ 3 ] ) + ".svg";

  // In out dirs.
  auto reader = JplaceReader();

  auto lhs = reader.read( from_file( argv[ 1 ] ) );
  auto rhs = reader.read( from_file( argv[ 2 ] ) );

  if( not compatible_trees( lhs, rhs ) ) {
    throw std::runtime_error{ "Trees are not compatible!" };
  }

  // get imbalance vectors for both samples
  auto imbalance_lhs = epca_imbalance_vector( lhs );
  auto imbalance_rhs = epca_imbalance_vector( rhs );
  assert( imbalance_lhs.size() == imbalance_rhs.size() );

  // calculate the difference
  for( size_t i = 0; i < imbalance_lhs.size(); ++i ) {
    imbalance_lhs[ i ] = ( imbalance_lhs[ i ] - imbalance_rhs[ i ] );
  }

  // map onto tree
  auto layout         = LayoutParameters{};
  layout.stroke.width = 15.0;
  layout.ladderize    = false;
  auto cmap           = ColorMap( color_list_spectral() );
  auto cnorm          = ColorNormalizationDiverging();
  cnorm.autoscale( imbalance_lhs );
  cnorm.make_centric();
  auto const colors = cmap( cnorm, imbalance_lhs );

  write_color_tree_to_svg_file( lhs.tree(), layout, colors, cmap, cnorm, outfile );

  // print_tree(lhs, "lhs.svg");
  // print_tree(rhs, "rhs.svg");

  LOG_INFO << "Finished";
  return 0;
}

#include "genesis/genesis.hpp"

genesis::sequence::SequenceSet read_any_seqfile( std::string const& file, bool const from_stdin=false )
{
  genesis::sequence::SequenceSet out_set;

  try {
    genesis::sequence::FastaReader().read( from_stdin
      ? genesis::utils::from_stream( std::cin )
      : genesis::utils::from_file( file ), out_set );
  } catch( std::exception& e ) {
    auto reader = genesis::sequence::PhylipReader();
    try {
      reader.read( from_stdin
        ? genesis::utils::from_stream( std::cin )
        : genesis::utils::from_file( file ), out_set );
    } catch( std::exception& e ) {
      out_set.clear();
      reader.mode( genesis::sequence::PhylipReader::Mode::kInterleaved );
      try {
        reader.read( from_stdin
          ? genesis::utils::from_stream( std::cin )
          : genesis::utils::from_file( file ), out_set );
      } catch( std::exception& e ) {
        throw std::invalid_argument {
          "Cannot parse sequence file(s): Invalid file format? (only phylip and fasta allowed)"
        };
      }
    }
  }

  return out_set;
}
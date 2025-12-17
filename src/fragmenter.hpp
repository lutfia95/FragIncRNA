#pragma once

#include "config.hpp"

#include <filesystem>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>

class Fragmenter
{
public:
    explicit Fragmenter(Config const & cfg);

    // Reads one reference (FASTA / FASTA.gz) and returns overlapping fragments.
    std::vector<seqan3::dna5_vector>
    fragment_reference(std::filesystem::path const & ref_path,
                       std::string const & ref_id) const;

private:
    Config cfg_;
};

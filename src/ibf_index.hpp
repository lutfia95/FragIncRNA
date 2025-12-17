#pragma once

#include "config.hpp"

#include <memory>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

class IBFIndex
{
public:
    IBFIndex(std::string ref_name,
             std::vector<seqan3::dna5_vector> const & fragments,
             Config const & cfg);

    [[nodiscard]] std::string const & ref_name() const noexcept
    {
        return ref_name_;
    }

    [[nodiscard]] seqan3::interleaved_bloom_filter<> & ibf()
    {
        return *ibf_;
    }

    [[nodiscard]] seqan3::interleaved_bloom_filter<> const & ibf() const
    {
        return *ibf_;
    }

    [[nodiscard]] std::size_t bin_count() const noexcept;

private:
    std::string ref_name_;
    std::unique_ptr<seqan3::interleaved_bloom_filter<>> ibf_;
};

#include "ibf_index.hpp"

#include "logger.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

#include <cereal/archives/binary.hpp>

//#include <seqan3/search/kmer/ungapped.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

IBFIndex::IBFIndex(std::string ref_name,
                   std::vector<seqan3::dna5_vector> const & fragments,
                   Config const & cfg)
    : ref_name_{std::move(ref_name)}
{
    using seqan3::bin_count;
    using seqan3::bin_size;
    using seqan3::hash_function_count;

    if (fragments.empty())
        throw std::runtime_error{"No fragments supplied for reference '" + ref_name_ + "'."};

    std::size_t max_len = 0;
    for (auto const & f : fragments)
        max_len = std::max<std::size_t>(max_len, f.size());

    if (cfg.kmer_size == 0 || cfg.kmer_size > max_len)
        throw std::runtime_error{
            "kmer_size must be > 0 and <= maximum fragment length for reference '" + ref_name_ + "'."};

    if (cfg.fpr <= 0.0 || cfg.fpr >= 1.0)
        throw std::runtime_error{
            "FPR must be in (0,1). Got " + std::to_string(cfg.fpr) +
            " for reference '" + ref_name_ + "'."};

    if (cfg.hash_functions == 0)
        throw std::runtime_error{
            "hash_functions must be >= 1 for reference '" + ref_name_ + "'."};

    std::size_t approx_kmers_per_frag =
        max_len >= cfg.kmer_size ? max_len - cfg.kmer_size + 1 : 1;

    // Compute Bloom filter size per bin from FPR, number of elements (n),
    // and number of hash functions (k):
    // p = (1 - exp(-k n / m))^k  =>  m = - (k * n) / ln(1 - p^(1/k))
    double n = static_cast<double>(approx_kmers_per_frag);
    double k = static_cast<double>(cfg.hash_functions);
    double p = cfg.fpr;

    double p_root = std::pow(p, 1.0 / k);
    double denom = std::log(1.0 - p_root);
    if (!std::isfinite(denom) || denom == 0.0)
        denom = -1.0; // Fallback to avoid NaN / inf

    double m_real = -k * n / denom;
    std::size_t bin_bits = static_cast<std::size_t>(std::ceil(m_real));

    // Ensure a minimum size to avoid degenerate bins
    bin_bits = std::max<std::size_t>(bin_bits, 1024u);

    std::size_t bins = fragments.size();

    ibf_ = std::make_unique<seqan3::interleaved_bloom_filter<>>(
        bin_count{bins},
        bin_size{bin_bits},
        hash_function_count{cfg.hash_functions}
    );

    auto hash_view = seqan3::views::kmer_hash(
        seqan3::ungapped{static_cast<uint8_t>(cfg.kmer_size)});


    std::size_t bin_idx = 0;
    for (auto const & frag : fragments)
    {
        for (auto h : frag | hash_view)
        {
            ibf_->emplace(h, seqan3::bin_index{bin_idx});
        }
        ++bin_idx;
    }
    Logger::print_stdout("Built IBF for reference '" + ref_name_ + "'", true);
    Logger::info("Built IBF for reference '" + ref_name_ + "' (" +
                 std::to_string(bins) + " bins, " +
                 std::to_string(bin_bits) + " bits per bin, " +
                 std::to_string(cfg.hash_functions) + " hash functions, " +
                 "target FPR=" + std::to_string(cfg.fpr) + ").");

    if (cfg.store_ibf)
    {
        std::filesystem::create_directories(cfg.output_dir);
        auto out_path = cfg.output_dir / (ref_name_ + ".ibf");
        std::ofstream os(out_path, std::ios::binary);
        if (!os)
            throw std::runtime_error{"Failed to open IBF output file: " + out_path.string()};
        cereal::BinaryOutputArchive archive(os);
        archive(*ibf_);
        Logger::info("Stored IBF for reference '" + ref_name_ +
                     "' to file: " + out_path.string());
    }
}

std::size_t IBFIndex::bin_count() const noexcept
{
    return ibf_ ? ibf_->bin_count() : 0u;
}

#include "config.hpp"
#include "fragmenter.hpp"
#include "ibf_index.hpp"
#include "logger.hpp"
#include "query_processor.hpp"

#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>

namespace fs = std::filesystem;

static bool is_fasta_like(fs::path const & p)
{
    std::string s = p.string();
    auto ends_with = [&s](std::string_view suf)
    {
        if (s.size() < suf.size())
            return false;
        return std::string_view{s}.substr(s.size() - suf.size()) == suf;
    };

    return ends_with(".fa") || ends_with(".fna") || ends_with(".fasta") ||
           ends_with(".fa.gz") || ends_with(".fna.gz") || ends_with(".fasta.gz");
}

int main(int argc, char ** argv)
{
    Config cfg;

    try
    {
        seqan3::argument_parser parser{"lncrna_mers", argc, argv};

        parser.info.short_description = "Per-reference IBF-based k-mer matcher for lncRNAs using SeqAn3.";
        parser.info.version = "0.1";

        parser.add_option(cfg.ref_dir,
                          'r', "ref-dir",
                          "Directory containing reference FASTA / FASTA.gz files.",
                          seqan3::option_spec::required);

        parser.add_option(cfg.query_file,
                          'q', "query-file",
                          "Query sequences in FASTA format.",
                          seqan3::option_spec::required);

        parser.add_option(cfg.fragment_size,
                  'f', "fragment-size",
                  "Fragment size used to split each reference (> 3).",
                  seqan3::option_spec::standard,
                  seqan3::arithmetic_range_validator<std::size_t>{4, std::numeric_limits<std::size_t>::max()});

        parser.add_option(cfg.kmer_size,
                  'k', "kmer-size",
                  "k-mer size for IBF construction and querying (default 14).",
                  seqan3::option_spec::standard,
                  seqan3::arithmetic_range_validator<std::size_t>{1, 32});

        parser.add_option(cfg.hash_functions,
                  'H', "hash-functions",
                  "Number of hash functions used in IBFs (default 3).",
                  seqan3::option_spec::standard,
                  seqan3::arithmetic_range_validator<std::size_t>{1, 32});

        parser.add_option(cfg.fpr,
                  'p', "fpr",
                  "Target false positive rate per IBF bin (default 0.01).",
                  seqan3::option_spec::standard,
                  seqan3::arithmetic_range_validator<double>{1e-9, 0.5});

        parser.add_option(cfg.hit_threshold,
                  't', "hit-threshold",
                  "Absolute k-mer hit threshold per (query, reference).",
                  seqan3::option_spec::standard,
                  seqan3::arithmetic_range_validator<std::uint64_t>{0, std::numeric_limits<std::uint64_t>::max()});

        parser.add_option(cfg.output_dir,
                          'O', "output-dir",
                          "Directory for all outputs (results, logs, IBFs, fragments).");

        parser.add_option(cfg.output_file,
                          '\0', "output-file",
                          "Result output file name (TSV) in combined mode.");

        parser.add_option(cfg.log_file,
                          '\0', "log-file",
                          "Log file name.");

        parser.add_flag(cfg.store_fragments,
                        '\0', "store-fragments",
                        "Store per-reference fragments as FASTA.");

        parser.add_flag(cfg.store_ibf,
                        '\0', "store-ibf",
                        "Store per-reference IBF indices as binary files.");

        parser.add_flag(cfg.cleanup_ibf,
                        '\0', "cleanup-ibf",
                        "Delete per-reference IBF files after query processing.");

        parser.add_option(cfg.single_results_writer,
                          '\0', "single-results-writer",
                          "true = one combined TSV; false = one results_<ref>.tsv per reference.",
                          seqan3::option_spec::standard);

        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cerr << "Argument parser error: " << ext.what() << '\n';
        return 1;
    }

    try
    {
        fs::create_directories(cfg.output_dir);
        auto log_path = cfg.output_dir / cfg.log_file;
        Logger::init(log_path.string());
        Logger::info("Starting lncrna_mers.");

        if (!fs::exists(cfg.ref_dir) || !fs::is_directory(cfg.ref_dir))
            throw std::runtime_error("Reference directory does not exist or is not a directory: " +
                                     cfg.ref_dir.string());

        if (cfg.kmer_size > cfg.fragment_size)
            throw std::runtime_error("kmer-size must be <= fragment-size.");

        if (cfg.fpr <= 0.0 || cfg.fpr >= 1.0)
            throw std::runtime_error("FPR must be in (0,1).");

        // Collect reference files
        std::vector<fs::path> ref_files;
        for (auto const & entry : fs::directory_iterator{cfg.ref_dir})
        {
            if (!entry.is_regular_file())
                continue;
            if (!is_fasta_like(entry.path()))
                continue;
            ref_files.push_back(entry.path());
        }

        if (ref_files.empty())
            throw std::runtime_error("No reference FASTA / FASTA.gz files found in: " +
                                     cfg.ref_dir.string());

        Logger::info("Found " + std::to_string(ref_files.size()) +
                     " reference files to index.");

        // Pre-scan queries to get IDs and count (for combined mode)
        std::vector<std::string> query_ids;
        {
            seqan3::sequence_file_input query_in{cfg.query_file};
            for (auto & rec : query_in)
                query_ids.emplace_back(rec.id());
        }
        std::size_t total_queries = query_ids.size();
        if (total_queries == 0)
            throw std::runtime_error("No queries found in file: " + cfg.query_file.string());

        Logger::info("Total queries: " + std::to_string(total_queries));

        // build ref IDs
        std::vector<std::string> ref_ids;
        ref_ids.reserve(ref_files.size());
        for (auto const & ref_path : ref_files)
        {
            auto stem      = ref_path.stem().string();
            auto inner     = fs::path{stem}.stem().string();
            std::string id = inner.empty() ? stem : inner;
            ref_ids.push_back(id);
        }

        // For combined mode, allocate results matrix
        std::vector<std::vector<RefResult>> results;
        if (cfg.single_results_writer)
        {
            results.assign(total_queries,
                           std::vector<RefResult>(ref_files.size()));
        }

        Fragmenter fragmenter{cfg};

        // Process each reference one-by-one
        for (std::size_t r = 0; r < ref_files.size(); ++r)
        {
            auto const & ref_path = ref_files[r];
            auto const & ref_id   = ref_ids[r];

            Logger::info("Processing reference: " + ref_path.string() +
                         " (id='" + ref_id + "')");

            auto fragments = fragmenter.fragment_reference(ref_path, ref_id);

            IBFIndex ibf_idx{ref_id, fragments, cfg};

            QueryProcessor qp{cfg, ibf_idx, ref_id};

            if (cfg.single_results_writer)
            {
                // fill column r of results
                qp.run_fill_results_col(r, results);
            }
            else
            {
                // write per-IBF TSV immediately
                auto out_path = cfg.output_dir / ("results_" + ref_id + ".tsv");
                qp.run_write_per_ibf(out_path);
            }

            // Delete IBF on disk immediately after use, if requested
            if (cfg.cleanup_ibf && cfg.store_ibf)
            {
                auto ibf_path = cfg.output_dir / (ref_id + ".ibf");
                std::error_code ec;
                if (fs::exists(ibf_path, ec))
                {
                    fs::remove(ibf_path, ec);
                    if (!ec)
                        Logger::info("Removed IBF file: " + ibf_path.string());
                    else
                        Logger::warn("Failed to remove IBF file: " + ibf_path.string() +
                                     " (" + ec.message() + ")");
                }
            }
        }

        // Combined mode: write final TSV once from results matrix
        if (cfg.single_results_writer)
        {
            auto out_path = cfg.output_dir / cfg.output_file;
            std::ofstream out(out_path);
            if (!out)
                throw std::runtime_error("Failed to open result output file: " + out_path.string());

            // header
            out << "query";
            for (std::size_t r = 0; r < ref_ids.size(); ++r)
            {
                out << '\t' << ref_ids[r] << "_count"
                    << '\t' << ref_ids[r] << "_pass"
                    << '\t' << ref_ids[r] << "_pct";
            }
            out << '\n';

            // rows
            for (std::size_t q = 0; q < total_queries; ++q)
            {
                out << query_ids[q];
                for (std::size_t r = 0; r < ref_ids.size(); ++r)
                {
                    auto const & rr = results[q][r];
                    out << '\t' << rr.count
                        << '\t' << (rr.pass ? 1 : 0)
                        << '\t' << std::fixed << std::setprecision(4) << rr.pct;
                }
                out << '\n';
            }

            Logger::info("Combined results written to: " + out_path.string());
        }

        Logger::info("lncrna_mers finished successfully.");
    }
    catch (std::exception const & e)
    {
        Logger::error(std::string{"Fatal error: "} + e.what());
        std::cerr << "Fatal error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}

// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "option_collation.hpp"

#include <string>
#include <iostream>
#include <cctype>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <iterator>
#include <algorithm>
#include <functional>
#include <utility>
#include <thread>
#include <sstream>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/lexical_cast.hpp>

#include "utils/path_utils.hpp"
#include "utils/read_stats.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/string_utils.hpp"
#include "utils/repeat_finder.hpp"
#include "utils/append.hpp"
#include "utils/maths.hpp"
#include "basics/phred.hpp"
#include "basics/genomic_region.hpp"
#include "basics/aligned_read.hpp"
#include "basics/ploidy_map.hpp"
#include "basics/trio.hpp"
#include "basics/pedigree.hpp"
#include "readpipe/read_pipe_fwd.hpp"
#include "core/tools/coretools.hpp"
#include "core/callers/caller_builder.hpp"
#include "logging/logging.hpp"
#include "io/region/region_parser.hpp"
#include "io/pedigree/pedigree_reader.hpp"
#include "io/variant/vcf_reader.hpp"
#include "io/variant/vcf_writer.hpp"
#include "exceptions/user_error.hpp"
#include "exceptions/program_error.hpp"
#include "exceptions/system_error.hpp"
#include "exceptions/missing_file_error.hpp"
#include "core/models/haplotype_likelihood_model.hpp"

namespace octopus { namespace options {

boost::optional<fs::path> get_output_path(const OptionMap& options);

// unsigned are banned from the option map to prevent user input errors, but once the option
// map is passed they are all safe
unsigned as_unsigned(const std::string& option, const OptionMap& options)
{
    return static_cast<unsigned>(options.at(option).as<int>());
}

bool is_run_command(const OptionMap& options)
{
    return options.count("help") == 0 && options.count("version") == 0;
}

bool is_debug_mode(const OptionMap& options)
{
    return options.count("debug") == 1;
}

bool is_trace_mode(const OptionMap& options)
{
    return options.count("trace") == 1;
}

namespace {

class InvalidWorkingDirectory : public UserError
{
    std::string do_where() const override
    {
        return "get_working_directory";
    }

    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "The working directory you specified ";
        ss << path_;
        ss << " does not exist";
        return ss.str();
    }

    std::string do_help() const override
    {
        return "enter a valid working directory";
    }

    fs::path path_;
public:
    InvalidWorkingDirectory(fs::path p) : path_ {std::move(p)} {}
};

fs::path get_working_directory(const OptionMap& options)
{
    if (options.count("working-directory") == 1) {
        auto result = expand_user_path(options.at("working-directory").as<fs::path>());
        if (!fs::exists(result) && !fs::is_directory(result)) {
            throw InvalidWorkingDirectory {result};
        }
        return result;
    }
    return fs::current_path();
}

fs::path resolve_path(const fs::path& path, const OptionMap& options)
{
    return ::octopus::resolve_path(path, get_working_directory(options));
}

struct Line
{
    std::string line_data;
    operator std::string() const
    {
        return line_data;
    }
};

std::istream& operator>>(std::istream& is, Line& data)
{
    std::getline(is, data.line_data);
    if (!data.line_data.empty() && data.line_data.back() == '\r') {
        data.line_data.pop_back();
    }
    return is;
}

std::vector<fs::path> extract_paths_from_file(const fs::path& file_path, const OptionMap& options)
{
    std::ifstream file {file_path.string()};
    assert(file.good());
    std::vector<fs::path> result {};
    std::transform(std::istream_iterator<Line>(file), std::istream_iterator<Line>(),
                   std::back_inserter(result), [] (const auto& line) { return line.line_data; });
    result.erase(std::remove_if(std::begin(result), std::end(result),
                                [] (const auto& path) { return path.empty(); }),
                 std::end(result));
    return result;
}

auto resolve_paths(const std::vector<fs::path>& paths, const OptionMap& options)
{
    std::vector<fs::path> result {};
    result.reserve(paths.size());
    std::transform(std::cbegin(paths), std::cend(paths), std::back_inserter(result),
                   [&] (const auto& path) { return resolve_path(path, options); });
    return result;
}

auto resolve_paths(const std::vector<std::string>& path_strings, const OptionMap& options)
{
    std::vector<fs::path> paths {std::cbegin(path_strings), std::cend(path_strings)};
    return resolve_paths(paths, options);
}

bool is_file_readable(const fs::path& path)
{
    std::ifstream tmp {path.string()};
    return tmp.good();
}

bool is_file_writable(const fs::path& path)
{
    if (!fs::exists(path.parent_path())) {
        return false;
    }
    std::ofstream test {path.string()};
    const bool result {test.is_open()};
    fs::remove(path);
    return result;
}

} // namespace

bool is_threading_allowed(const OptionMap& options)
{
    unsigned num_threads {1};
    if (options.count("threads") == 1) {
        num_threads = as_unsigned("threads", options);
    }
    return num_threads != 1;
}

boost::optional<unsigned> get_num_threads(const OptionMap& options)
{
    unsigned num_threads {1};
    if (options.count("threads") == 1) {
        num_threads = as_unsigned("threads", options);
    }
    if (num_threads > 0) return num_threads;
    return boost::none;
}

ExecutionPolicy get_thread_execution_policy(const OptionMap& options)
{
    if (options.count("threads") == 1) {
        if (options.at("threads").as<int>() == 0) {
            return ExecutionPolicy::par;
        } else {
            return ExecutionPolicy::seq;
        }
    } else {
        return ExecutionPolicy::seq;
    }
}

MemoryFootprint get_target_read_buffer_size(const OptionMap& options)
{
    return options.at("target-read-buffer-footprint").as<MemoryFootprint>();
}

boost::optional<fs::path> get_debug_log_file_name(const OptionMap& options)
{
    if (options.count("debug") == 1) {
        return resolve_path(options.at("debug").as<fs::path>(), options);
    } else {
        return boost::none;
    }
}

boost::optional<fs::path> get_trace_log_file_name(const OptionMap& options)
{
    if (options.count("trace") == 1) {
        return resolve_path(options.at("trace").as<fs::path>(), options);
    } else {
        return boost::none;
    }
}

bool is_fast_mode(const OptionMap& options)
{
    return options.at("fast").as<bool>() || options.at("very-fast").as<bool>();
}

bool is_very_fast_mode(const OptionMap& options)
{
    return options.at("very-fast").as<bool>();
}

ReferenceGenome make_reference(const OptionMap& options)
{
    const fs::path input_path {options.at("reference").as<std::string>()};
    auto resolved_path = resolve_path(input_path, options);
    const auto ref_cache_size = options.at("max-reference-cache-footprint").as<MemoryFootprint>().num_bytes();
    try {
        return octopus::make_reference(std::move(resolved_path),
                                       ref_cache_size,
                                       is_threading_allowed(options));
    } catch (MissingFileError& e) {
        e.set_location_specified("the command line option --reference");
        throw;
    } catch (...) {
        throw;
    }
}

InputRegionMap make_search_regions(const std::vector<GenomicRegion>& regions)
{
    std::map<ContigName, std::deque<GenomicRegion>> contig_mapped_regions {};
    for (const auto& region : regions) {
        contig_mapped_regions[region.contig_name()].push_back(region);
    }
    InputRegionMap result {};
    result.reserve(contig_mapped_regions.size());
    for (auto& p : contig_mapped_regions) {
        std::sort(std::begin(p.second), std::end(p.second));
        auto covered_contig_regions = extract_covered_regions(p.second);
        result.emplace(std::piecewise_construct,
                       std::forward_as_tuple(p.first),
                       std::forward_as_tuple(std::make_move_iterator(std::begin(covered_contig_regions)),
                                             std::make_move_iterator(std::end(covered_contig_regions))));
    }
    return result;
}

InputRegionMap extract_search_regions(const ReferenceGenome& reference)
{
    return make_search_regions(get_all_contig_regions(reference));
}

auto get_unskipped(const MappableFlatSet<GenomicRegion>& regions, const MappableFlatSet<GenomicRegion>& skips)
{
    if (skips.empty()) return regions;
    MappableFlatSet<GenomicRegion> result {};
    for (const auto& region : regions) {
        const auto overlapped = skips.overlap_range(region);
        if (empty(overlapped)) {
            result.emplace(region);
        } else if (!contains(overlapped.front(), region)) {
            if (begins_before(region, overlapped.front())) {
                result.emplace(left_overhang_region(region, overlapped.front()));
            }
            auto intervening_chunks = extract_intervening_regions(overlapped);
            using std::make_move_iterator;
            result.insert(make_move_iterator(std::begin(intervening_chunks)),
                          make_move_iterator(std::end(intervening_chunks)));
            if (ends_before(overlapped.back(), region)) {
                result.emplace(right_overhang_region(region, overlapped.back()));
            }
        }
    }
    result.shrink_to_fit();
    return result;
}

InputRegionMap extract_search_regions(const std::vector<GenomicRegion>& regions,
                                      std::vector<GenomicRegion>& skip_regions)
{
    auto input_regions = make_search_regions(regions);
    const auto skipped = make_search_regions(skip_regions);
    InputRegionMap result {input_regions.size()};
    for (auto& p : input_regions) {
        if (skipped.count(p.first) == 1) {
            result.emplace(p.first, get_unskipped(std::move(p.second), skipped.at(p.first)));
        } else {
            result.emplace(p.first, std::move(p.second));
        }
    }
    for (auto it = std::begin(result); it != std::end(result); ) {
        if (it->second.empty()) {
            it = result.erase(it);
        } else {
            ++it;
        }
    }
    for (auto& p : result) {
        p.second.shrink_to_fit();
    }
    return result;
}

InputRegionMap extract_search_regions(const ReferenceGenome& reference,
                                      std::vector<GenomicRegion>& skip_regions)
{
    return extract_search_regions(get_all_contig_regions(reference), skip_regions);
}
        
std::vector<GenomicRegion> parse_regions(const std::vector<std::string>& unparsed_regions,
                                         const ReferenceGenome& reference)
{
    std::vector<GenomicRegion> result {};
    result.reserve(unparsed_regions.size());
    for (const auto& unparsed_region : unparsed_regions) {
        result.push_back(io::parse_region(unparsed_region, reference));
    }
    return result;
}

auto transform_to_zero_based(std::vector<GenomicRegion>&& one_based_regions)
{
    std::vector<GenomicRegion> result {};
    result.reserve(one_based_regions.size());
    for (auto&& region : one_based_regions) {
        if (region.begin() > 0) {
            result.push_back(shift(std::move(region), -1));
        } else {
            result.push_back(std::move(region));
        }
    }
    return result;
}

auto transform_to_zero_based(InputRegionMap::mapped_type&& one_based_regions)
{
    MappableFlatSet<GenomicRegion> result {};
    for (auto&& region : one_based_regions) {
        result.insert(shift(std::move(region), -1));
    }
    return result;
}

auto transform_to_zero_based(InputRegionMap&& one_based_search_regions)
{
    InputRegionMap result {one_based_search_regions.size()};
    for (auto& p : one_based_search_regions) {
        result.emplace(p.first, transform_to_zero_based(std::move(p.second)));
    }
    return result;
}

class MissingRegionPathFile : public MissingFileError
{
    std::string do_where() const override
    {
        return "get_search_regions";
    }
public:
    MissingRegionPathFile(fs::path p) : MissingFileError {std::move(p), "region path"} {};
};

InputRegionMap get_search_regions(const OptionMap& options, const ReferenceGenome& reference)
{
    using namespace utils;
    std::vector<GenomicRegion> skip_regions {};
    if (options.count("skip-regions") == 1) {
        const auto& region_strings = options.at("skip-regions").as<std::vector<std::string>>();
        append(parse_regions(region_strings, reference), skip_regions);
    }
    if (options.count("skip-regions-file") == 1) {
        const auto& input_path = options.at("skip-regions-file").as<fs::path>();
        auto resolved_path = resolve_path(input_path, options);
        if (!fs::exists(resolved_path)) {
            MissingRegionPathFile e {resolved_path};
            e.set_location_specified("the command line option '--skip-regions-file'");
            throw e;
        }
        auto regions = io::extract_regions(resolved_path, reference, io::NonreferenceContigPolicy::ignore);
        if (regions.empty()) {
            logging::WarningLogger log {};
            stream(log) << "The regions path file you specified " << resolved_path
                << " in the command line option '--skip-regions-file' is empty";
        }
        append(std::move(regions), skip_regions);
    }
    if (options.at("one-based-indexing").as<bool>()) {
        skip_regions = transform_to_zero_based(std::move(skip_regions));
    }
    if (options.count("regions") == 0 && options.count("regions-file") == 0) {
        if (options.count("regenotype") == 1) {
            // TODO: only extract regions in the regenotype VCF
        }
        return extract_search_regions(reference, skip_regions);
    }
    std::vector<GenomicRegion> input_regions {};
    if (options.count("regions") == 1) {
        const auto& region_strings = options.at("regions").as<std::vector<std::string>>();
        append(parse_regions(region_strings, reference), input_regions);
    }
    if (options.count("regions-file") == 1) {
        const auto& input_path = options.at("regions-file").as<fs::path>();
        auto resolved_path = resolve_path(input_path, options);
        if (!fs::exists(resolved_path)) {
            MissingRegionPathFile e {resolved_path};
            e.set_location_specified("the command line option '--regions-file'");
            throw e;
        }
        auto regions = io::extract_regions(resolved_path, reference);
        if (regions.empty()) {
            logging::WarningLogger log {};
            stream(log) << "The regions path file you specified " << resolved_path
                << " in the command line option '--skip-regions-file' is empty";
        }
        append(std::move(regions), input_regions);
    }
    auto result = extract_search_regions(input_regions, skip_regions);
    if (options.at("one-based-indexing").as<bool>()) {
        return transform_to_zero_based(std::move(result));
    }
    return result;
}

ContigOutputOrder get_contig_output_order(const OptionMap& options)
{
    return options.at("contig-output-order").as<ContigOutputOrder>();
}

bool ignore_unmapped_contigs(const OptionMap& options)
{
    return options.at("ignore-unmapped-contigs").as<bool>();
}

boost::optional<std::vector<SampleName>> get_user_samples(const OptionMap& options)
{
    if (options.count("samples") == 1) {
        return options.at("samples").as<std::vector<SampleName>>();
    }
    return boost::none;
}

class MissingReadPathFile : public MissingFileError
{
    std::string do_where() const override
    {
        return "get_read_paths";
    }
public:
    MissingReadPathFile(fs::path p) : MissingFileError {std::move(p), "read path"} {};
};

void log_and_remove_duplicates(std::vector<fs::path>& paths)
{
    std::sort(std::begin(paths), std::end(paths));
    const auto first_duplicate = std::adjacent_find(std::begin(paths), std::end(paths));
    if (first_duplicate != std::end(paths)) {
        std::deque<fs::path> duplicates {};
        for (auto duplicate_itr = first_duplicate; duplicate_itr != std::end(paths); ) {
            duplicates.push_back(*duplicate_itr);
            duplicate_itr = std::adjacent_find(std::find_if(std::next(duplicate_itr, 2), std::end(paths),
                                                            [=] (const auto& path) { return path != *duplicate_itr; }),
                                               std::end(paths));
        }
        const auto num_paths = paths.size();
        paths.erase(std::unique(first_duplicate, std::end(paths)), std::end(paths));
        const auto num_unique_paths = paths.size();
        const auto num_duplicate_paths = num_paths - num_unique_paths;
        logging::WarningLogger warn_log {};
        auto warn_log_stream = stream(warn_log);
        warn_log_stream << "Ignoring " << num_duplicate_paths << " duplicate read path";
        if (num_duplicate_paths > 1) {
            warn_log_stream << 's';
        }
        warn_log_stream << ": ";
        std::for_each(std::cbegin(duplicates), std::prev(std::cend(duplicates)), [&] (const auto& path) {
            warn_log_stream << path << ", ";
        });
        warn_log_stream << duplicates.back();
        if (num_duplicate_paths > duplicates.size()) {
            warn_log_stream << " (showing unique duplicates)";
        }
    }
}

std::vector<fs::path> get_read_paths(const OptionMap& options)
{
    using namespace utils;
    std::vector<fs::path> result {};
    if (options.count("reads") == 1) {
        auto resolved_paths = resolve_paths(options.at("reads").as<std::vector<fs::path>>(), options);
        append(std::move(resolved_paths), result);
    }
    if (options.count("reads-file") == 1) {
        const fs::path input_path {options.at("reads-file").as<fs::path>()};
        auto resolved_path = resolve_path(input_path, options);
        if (!fs::exists(resolved_path)) {
            MissingReadPathFile e {resolved_path};
            e.set_location_specified("the command line option '--reads-file'");
            throw e;
        }
        auto paths = extract_paths_from_file(resolved_path, options);
        auto resolved_paths = resolve_paths(paths, options);
        if (resolved_paths.empty()) {
            logging::WarningLogger log {};
            stream(log) << "The read path file you specified " << resolved_path
                        << " in the command line option '--reads-file' is empty";
        }
        append(std::move(resolved_paths), result);
    }
    log_and_remove_duplicates(result);
    return result;
}

ReadManager make_read_manager(const OptionMap& options)
{
    auto read_paths = get_read_paths(options);
    const auto max_open_files = as_unsigned("max-open-read-files", options);
    return ReadManager {std::move(read_paths), max_open_files};
}

bool allow_assembler_generation(const OptionMap& options)
{
    return options.at("assembly-candidate-generator").as<bool>() && !is_fast_mode(options);
}

auto make_read_transformers(const OptionMap& options)
{
    using namespace octopus::readpipe;
    ReadTransformer prefilter_transformer {}, postfilter_transformer {};
    prefilter_transformer.add(CapitaliseBases {});
    prefilter_transformer.add(CapBaseQualities {125});
    if (options.at("read-transforms").as<bool>()) {
        if (options.count("mask-low-quality-tails") == 1) {
            const auto threshold = static_cast<AlignedRead::BaseQuality>(as_unsigned("mask-low-quality-tails", options));
            prefilter_transformer.add(MaskLowQualityTails {threshold});
        }
        if (options.at("soft-clip-masking").as<bool>()) {
            const auto boundary_size = as_unsigned("mask-soft-clipped-boundary-bases", options);
            if (boundary_size > 0) {
                if (options.count("soft-clip-mask-threshold") == 1) {
                    const auto threshold = static_cast<AlignedRead::BaseQuality>(as_unsigned("soft-clip-mask-threshold", options));
                    prefilter_transformer.add(MaskLowQualitySoftClippedBoundaryBases {boundary_size, threshold});
                } else if (allow_assembler_generation(options)) {
                    prefilter_transformer.add(MaskLowQualitySoftClippedBoundaryBases {boundary_size, 3});
                    prefilter_transformer.add(MaskLowAverageQualitySoftClippedTails {10, 5});
                    prefilter_transformer.add(MaskClippedDuplicatedBases {});
                } else {
                    prefilter_transformer.add(MaskSoftClippedBoundraryBases {boundary_size});
                }
            } else {
                if (options.count("soft-clip-mask-threshold") == 1) {
                    const auto threshold = static_cast<AlignedRead::BaseQuality>(as_unsigned("soft-clip-mask-threshold", options));
                    prefilter_transformer.add(MaskLowQualitySoftClippedBases {threshold});
                } else if (allow_assembler_generation(options)) {
                    prefilter_transformer.add(MaskLowQualitySoftClippedBases {3});
                    prefilter_transformer.add(MaskLowAverageQualitySoftClippedTails {10, 5});
                    prefilter_transformer.add(MaskClippedDuplicatedBases {});
                } else {
                    prefilter_transformer.add(MaskSoftClipped {});
                }
            }
        }
        if (options.at("adapter-masking").as<bool>()) {
            prefilter_transformer.add(MaskAdapters {});
            postfilter_transformer.add(MaskTemplateAdapters {});
        }
        if (options.at("overlap-masking").as<bool>()) {
            postfilter_transformer.add(MaskStrandOfDuplicatedBases {});
        }
        prefilter_transformer.shrink_to_fit();
        postfilter_transformer.shrink_to_fit();
    }
    return std::make_pair(std::move(prefilter_transformer), std::move(postfilter_transformer));
}

bool is_read_filtering_enabled(const OptionMap& options)
{
    return options.at("read-filtering").as<bool>();
}

auto make_read_filterer(const OptionMap& options)
{
    using std::make_unique;
    using namespace octopus::readpipe;
    using ReadFilterer = ReadPipe::ReadFilterer;
    
    ReadFilterer result {};
    
    // these filters are mandatory
    result.add(make_unique<HasValidBaseQualities>());
    result.add(make_unique<HasWellFormedCigar>());
    
    if (!is_read_filtering_enabled(options)) {
        return result;
    }
    if (!options.at("consider-unmapped-reads").as<bool>()) {
        result.add(make_unique<IsMapped>());
    }
    
    const auto min_mapping_quality = as_unsigned("min-mapping-quality", options);
    const auto min_base_quality    = as_unsigned("good-base-quality", options);
    const auto min_good_bases      = as_unsigned("min-good-bases", options);
    
    if (min_mapping_quality > 0) {
        result.add(make_unique<IsGoodMappingQuality>(min_mapping_quality));
    }
    if (min_base_quality > 0 && min_good_bases > 0) {
        result.add(make_unique<HasSufficientGoodQualityBases>(min_base_quality, min_good_bases));
    }
    if (min_base_quality > 0 && options.count("min-good-base-fraction") == 1) {
        auto min_good_base_fraction = options.at("min-good-base-fraction").as<double>();
        result.add(make_unique<HasSufficientGoodBaseFraction>(min_base_quality, min_good_base_fraction));
    }
    if (options.count("min-read-length") == 1) {
        result.add(make_unique<IsShort>(as_unsigned("min-read-length", options)));
    }
    if (options.count("max-read-length") == 1) {
        result.add(make_unique<IsLong>(as_unsigned("max-read-length", options)));
    }
    if (!options.at("allow-marked-duplicates").as<bool>()) {
        result.add(make_unique<IsNotMarkedDuplicate>());
    }
    if (!options.at("allow-octopus-duplicates").as<bool>()) {
        result.add(make_unique<IsNotDuplicate<ReadFilterer::ReadIterator>>());
    }
    if (!options.at("allow-qc-fails").as<bool>()) {
        result.add(make_unique<IsNotMarkedQcFail>());
    }
    if (options.at("no-secondary-alignments").as<bool>()) {
        result.add(make_unique<IsNotSecondaryAlignment>());
    }
    if (options.at("no-supplementary-alignments").as<bool>()) {
        result.add(make_unique<IsNotSupplementaryAlignment>());
    }
    if (!options.at("consider-reads-with-unmapped-segments").as<bool>()) {
        result.add(make_unique<IsNextSegmentMapped>());
        result.add(make_unique<IsProperTemplate>());
    }
    if (!options.at("consider-reads-with-distant-segments").as<bool>()) {
        result.add(make_unique<IsLocalTemplate>());
    }
    if (options.at("no-adapter-contaminated-reads").as<bool>()) {
        result.add(make_unique<IsNotContaminated>());
    }
    result.shrink_to_fit();
    return result;
}

bool is_downsampling_enabled(const OptionMap& options)
{
    return is_read_filtering_enabled(options) && !options.at("disable-downsampling").as<bool>();
}

boost::optional<readpipe::Downsampler> make_downsampler(const OptionMap& options)
{
    if (is_downsampling_enabled(options)) {
        using namespace octopus::readpipe;
        const auto max_coverage    = as_unsigned("downsample-above", options);
        const auto target_coverage = as_unsigned("downsample-target", options);
        return Downsampler {max_coverage, target_coverage};
    }
    return boost::none;
}

ReadPipe make_read_pipe(ReadManager& read_manager, std::vector<SampleName> samples, const OptionMap& options)
{
    auto transformers = make_read_transformers(options);
    if (transformers.second.num_transforms() > 0) {
        return ReadPipe {read_manager, std::move(transformers.first), make_read_filterer(options),
                         std::move(transformers.second), make_downsampler(options), std::move(samples)};
    } else {
        return ReadPipe {read_manager, std::move(transformers.first), make_read_filterer(options),
                         make_downsampler(options), std::move(samples)};
    }
}

auto get_default_inclusion_predicate()
{
    return coretools::DefaultInclusionPredicate {};
}

auto get_default_inclusion_predicate(const OptionMap& options) noexcept
{
    using namespace coretools;
    using InclusionPredicate = DynamicCigarScanner::Options::InclusionPredicate;
    const auto caller = options.at("caller").as<std::string>();
    if (caller == "cancer") {
        // TODO: specialise for this case; we need to be careful about low frequency somatics.
        return InclusionPredicate {get_default_inclusion_predicate()};
    } else {
        return InclusionPredicate {get_default_inclusion_predicate()};
    }
}

auto get_default_match_predicate() noexcept
{
    return coretools::DefaultMatchPredicate {};
}

class MissingSourceVariantFile : public MissingFileError
{
    std::string do_where() const override
    {
        return "make_variant_generator_builder";
    }
public:
    MissingSourceVariantFile(fs::path p) : MissingFileError {std::move(p), "source variant"} {};
};

class ConflictingSourceVariantFile : public UserError
{
    std::string do_where() const override
    {
        return "make_variant_generator_builder";
    }
    
    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "The source variant file you specified " << source_;
        ss << " conflicts with the output file " << output_;
        return ss.str();
    }
    
    std::string do_help() const override
    {
        return "Specify a unique output file";
    }
    
    fs::path source_, output_;
public:
    ConflictingSourceVariantFile(fs::path source, fs::path output)
    : source_ {std::move(source)}
    , output_ {std::move(output)}
    {}
};

struct DefaultRepeatGenerator
{
    std::vector<GenomicRegion> operator()(const ReferenceGenome& reference, GenomicRegion region) const
    {
        return find_repeat_regions(reference, region);
    }
};

auto get_max_expected_heterozygosity(const OptionMap& options)
{
    const auto snp_heterozygosity = options.at("snp-heterozygosity").as<float>();
    const auto indel_heterozygosity = options.at("indel-heterozygosity").as<float>();
    const auto heterozygosity = snp_heterozygosity + indel_heterozygosity;
    const auto heterozygosity_stdev = options.at("snp-heterozygosity-stdev").as<float>();
    return std::min(static_cast<double>(heterozygosity + 2 * heterozygosity_stdev), 0.9999);
}

auto make_variant_generator_builder(const OptionMap& options)
{
    using namespace coretools;
    
    logging::WarningLogger warning_log {};
    logging::ErrorLogger log {};
    
    VariantGeneratorBuilder result {};
    const bool use_assembler {allow_assembler_generation(options)};
    
    if (options.at("raw-cigar-candidate-generator").as<bool>()) {
        if (options.count("min-supporting-reads") == 1) {
            CigarScanner::Options scanner_options {};
            scanner_options.min_base_quality = as_unsigned("min-base-quality", options);
            scanner_options.min_support = as_unsigned("min-supporting-reads", options);
            if (scanner_options.min_support == 0) {
                warning_log << "The option --min_supporting_reads was set to 0 - assuming this is a typo and setting to 1";
                ++scanner_options.min_support;
            }
            result.set_cigar_scanner(scanner_options);
        } else {
            DynamicCigarScanner::Options scanner_options {};
            scanner_options.include = get_default_inclusion_predicate(options);
            scanner_options.match = get_default_match_predicate();
            scanner_options.use_clipped_coverage_tracking = true;
            scanner_options.repeat_region_generator = DefaultRepeatGenerator {};
            DynamicCigarScanner::Options::MisalignmentParameters misalign_params {};
            misalign_params.max_expected_mutation_rate = get_max_expected_heterozygosity(options);
            misalign_params.snv_threshold = as_unsigned("min-base-quality", options);
            if (use_assembler) {
                misalign_params.indel_penalty = 1.5;
                misalign_params.clip_penalty = 2;
                misalign_params.min_ln_prob_correctly_aligned = std::log(0.005);
            }
            scanner_options.misalignment_parameters = misalign_params;
            result.set_dynamic_cigar_scanner(std::move(scanner_options));
        }
    }
    if (use_assembler) {
        LocalReassembler::Options reassembler_options {};
        const auto kmer_sizes = options.at("kmer-sizes").as<std::vector<int>>();
        reassembler_options.kmer_sizes.assign(std::cbegin(kmer_sizes), std::cend(kmer_sizes));
        if (options.count("assembler-mask-base-quality") == 1) {
            reassembler_options.mask_threshold = as_unsigned("assembler-mask-base-quality", options);
        }
        reassembler_options.execution_policy = get_thread_execution_policy(options);
        reassembler_options.num_fallbacks = as_unsigned("num-fallback-kmers", options);
        reassembler_options.fallback_interval_size = as_unsigned("fallback-kmer-gap", options);
        reassembler_options.bin_size = as_unsigned("max-region-to-assemble", options);
        reassembler_options.bin_overlap = as_unsigned("max-assemble-region-overlap", options);
        reassembler_options.min_kmer_observations = as_unsigned("min-kmer-prune", options);
        reassembler_options.max_bubbles = as_unsigned("max-bubbles", options);
        reassembler_options.min_bubble_score = options.at("min-bubble-score").as<double>();
        reassembler_options.max_variant_size = as_unsigned("max-variant-size", options);
        result.set_local_reassembler(std::move(reassembler_options));
    }
    if (options.count("source-candidates") == 1) {
        const auto output_path = get_output_path(options);
        const auto input_paths = options.at("source-candidates").as<std::vector<fs::path>>();
        for (const auto& input_path : input_paths) {
            auto resolved_source_path = resolve_path(input_path, options);
            if (!fs::exists(resolved_source_path)) {
                throw MissingSourceVariantFile {input_path};
            }
            if (output_path && resolved_source_path == *output_path) {
                throw ConflictingSourceVariantFile {std::move(resolved_source_path), *output_path};
            }
            result.add_vcf_extractor(std::move(resolved_source_path));
        }
    }
    if (options.count("regenotype") == 1) {
        auto regenotype_path = options.at("regenotype").as<fs::path>();
        if (options.count("source-candidates") == 1) {
            fs::path input_path {options.at("source-candidates").as<std::string>()};
            if (regenotype_path != input_path) {
                warning_log << "Running in regenotype mode but given a different source variant file";
            }
            return result;
        }
        auto resolved_regenotype_path = resolve_path(regenotype_path, options);
        if (!fs::exists(resolved_regenotype_path)) {
            throw MissingSourceVariantFile {resolved_regenotype_path};
        }
        const auto output_path = get_output_path(options);
        if (output_path && resolved_regenotype_path == *output_path) {
            throw ConflictingSourceVariantFile {std::move(resolved_regenotype_path), *output_path};
        }
        result.add_vcf_extractor(std::move(resolved_regenotype_path));
    }
    
    return result;
}

struct ContigPloidyLess
{
    bool operator()(const ContigPloidy& lhs, const ContigPloidy& rhs) const noexcept
    {
        if (lhs.sample) {
            if (rhs.sample && *lhs.sample != *rhs.sample) {
                return *lhs.sample < *rhs.sample;
            } else {
                return true;
            }
        } else if (rhs.sample) {
            return false;
        }
        return lhs.contig == rhs.contig ? lhs.ploidy < rhs.ploidy : lhs.contig < rhs.contig;
    }
};

struct ContigPloidyEqual
{
    bool operator()(const ContigPloidy& lhs, const ContigPloidy& rhs) const noexcept
    {
        return lhs.sample == rhs.sample && lhs.contig == rhs.contig && lhs.ploidy == rhs.ploidy;
    }
};

struct ContigPloidyAmbiguous
{
    bool operator()(const ContigPloidy& lhs, const ContigPloidy& rhs) const noexcept
    {
        if (lhs.sample && rhs.sample) {
            return *lhs.sample == *rhs.sample && lhs.contig == rhs.contig;
        } else if (!(lhs.sample || rhs.sample)) {
            return lhs.contig == rhs.contig;
        }
        return false;
    }
};

class AmbiguousPloidy : public UserError
{
    std::string do_where() const override
    {
        return "make_caller_factory";
    }
    
    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "The are contigs with ambiguous ploidy: ";
        for (auto it = std::cbegin(ploidies_), end = std::cend(ploidies_); it != end;) {
            it = std::adjacent_find(it, std::cend(ploidies_), ContigPloidyAmbiguous {});
            if (it != std::cend(ploidies_)) {
                const auto it2 = std::find_if(std::next(it), std::cend(ploidies_),
                                              [=] (const auto& cp) {
                                                  return ContigPloidyAmbiguous{}(*it, cp);
                                              });
                std::ostringstream ss {};
                std::copy(it, it2, std::ostream_iterator<ContigPloidy> {ss, " "});
                it = it2;
            }
        }
        return ss.str();
    }
    
    std::string do_help() const override
    {
        return "Ensure ploidies are specified only once per sample or per sample contig";
    }
    
    std::vector<ContigPloidy> ploidies_;

public:
    AmbiguousPloidy(std::vector<ContigPloidy> ploidies) : ploidies_ {ploidies} {}
};

void remove_duplicate_ploidies(std::vector<ContigPloidy>& contig_ploidies)
{
    std::sort(std::begin(contig_ploidies), std::end(contig_ploidies), ContigPloidyLess {});
    auto itr = std::unique(std::begin(contig_ploidies), std::end(contig_ploidies), ContigPloidyEqual {});
    contig_ploidies.erase(itr, std::end(contig_ploidies));
}

bool has_ambiguous_ploidies(const std::vector<ContigPloidy>& contig_ploidies)
{
    return std::adjacent_find(std::cbegin(contig_ploidies), std::cend(contig_ploidies),
                              ContigPloidyAmbiguous {}) != std::cend(contig_ploidies);
}

class MissingPloidyFile : public MissingFileError
{
    std::string do_where() const override
    {
        return "get_ploidy_map";
    }
public:
    MissingPloidyFile(fs::path p) : MissingFileError {std::move(p), "ploidy"} {};
};

PloidyMap get_ploidy_map(const OptionMap& options)
{
    std::vector<ContigPloidy> flat_plodies {};
    if (options.count("contig-ploidies-file") == 1) {
        const fs::path input_path {options.at("contig-ploidies-file").as<std::string>()};
        const auto resolved_path = resolve_path(input_path, options);
        if (!fs::exists(resolved_path)) {
            throw MissingPloidyFile {input_path};
        }
        std::ifstream file {resolved_path.string()};
        std::transform(std::istream_iterator<Line>(file), std::istream_iterator<Line>(),
                       std::back_inserter(flat_plodies), [] (const Line& line) {
            std::istringstream ss {line.line_data};
            ContigPloidy result {};
            ss >> result;
            return result;
        });
    }
    if (options.count("contig-ploidies") == 1) {
        utils::append(options.at("contig-ploidies").as<std::vector<ContigPloidy>>(), flat_plodies);
    }
    remove_duplicate_ploidies(flat_plodies);
    if (has_ambiguous_ploidies(flat_plodies)) {
        throw AmbiguousPloidy {flat_plodies};
    }
    PloidyMap result {as_unsigned("organism-ploidy", options)};
    for (const auto& p : flat_plodies) {
        if (p.sample) {
            result.set(*p.sample, p.contig, p.ploidy);
        } else {
            result.set(p.contig, p.ploidy);
        }
    }
    return result;
}

bool call_sites_only(const OptionMap& options)
{
    return options.at("sites-only").as<bool>();
}

auto get_extension_policy(const OptionMap& options)
{
    using ExtensionPolicy = HaplotypeGenerator::Builder::Policies::Extension;
    switch (options.at("extension-level").as<ExtensionLevel>()) {
        case ExtensionLevel::conservative: return ExtensionPolicy::conservative;
        case ExtensionLevel::normal: return ExtensionPolicy::normal;
        case ExtensionLevel::optimistic: return ExtensionPolicy::optimistic;
        case ExtensionLevel::aggressive: return ExtensionPolicy::aggressive;
        default: return ExtensionPolicy::normal; // to stop GCC warning
    }
}

auto get_lagging_policy(const OptionMap& options)
{
    using LaggingPolicy = HaplotypeGenerator::Builder::Policies::Lagging;
    if (is_fast_mode(options)) return LaggingPolicy::none;
    switch (options.at("phasing-level").as<PhasingLevel>()) {
        case PhasingLevel::conservative: return LaggingPolicy::conservative;
        case PhasingLevel::moderate: return LaggingPolicy::moderate;
        case PhasingLevel::normal: return LaggingPolicy::normal;
        case PhasingLevel::aggressive: return LaggingPolicy::aggressive;
        default: return LaggingPolicy::none;
    }
}

auto get_max_haplotypes(const OptionMap& options)
{
    if (is_fast_mode(options)) {
        return 50u;
    } else {
        return as_unsigned("max-haplotypes", options);
    }
}

auto get_max_expected_log_allele_count_per_base(const OptionMap& options)
{
    const auto snp_heterozygosity = options.at("snp-heterozygosity").as<float>();
    const auto indel_heterozygosity = options.at("indel-heterozygosity").as<float>();
    const auto heterozygosity = snp_heterozygosity + indel_heterozygosity;
    const auto snp_heterozygosity_stdev = options.at("snp-heterozygosity-stdev").as<float>();
    const auto max_log_allele_count_per_base = heterozygosity + 8 * snp_heterozygosity_stdev;
    return max_log_allele_count_per_base;
}

auto get_max_indicator_join_distance() noexcept
{
    return HaplotypeLikelihoodModel{}.pad_requirement();
}

auto get_min_flank_pad() noexcept
{
    return 2 * (2 * HaplotypeLikelihoodModel{}.pad_requirement() - 1);
}

auto make_haplotype_generator_builder(const OptionMap& options)
{
    const auto lagging_policy    = get_lagging_policy(options);
    const auto max_haplotypes    = get_max_haplotypes(options);
    const auto holdout_limit     = as_unsigned("haplotype-holdout-threshold", options);
    const auto overflow_limit    = as_unsigned("haplotype-overflow", options);
    const auto max_holdout_depth = as_unsigned("max-holdout-depth", options);
    return HaplotypeGenerator::Builder().set_extension_policy(get_extension_policy(options))
    .set_target_limit(max_haplotypes).set_holdout_limit(holdout_limit).set_overflow_limit(overflow_limit)
    .set_lagging_policy(lagging_policy).set_max_holdout_depth(max_holdout_depth)
    .set_max_indicator_join_distance(get_max_indicator_join_distance())
    .set_max_expected_log_allele_count_per_base(get_max_expected_log_allele_count_per_base(options))
    .set_min_flank_pad(get_min_flank_pad());
}

boost::optional<Pedigree> get_pedigree(const OptionMap& options)
{
	if (options.count("pedigree") == 1) {
		const auto ped_file = resolve_path(options.at("pedigree").as<fs::path>(), options);
		return io::read_pedigree(ped_file);
	} else {
		return boost::none;
	}
}

class BadTrioSampleSet : public UserError
{
    std::string do_where() const override
    {
        return "make_trio";
    }
    
    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "Trio calling requires exactly 3 samples but "
           << num_samples_
           << " where provided";
        return ss.str();
    }
    
    std::string do_help() const override
    {
        return "Ensure only three samples are present; if the read files contain more than"
                " this then explicitly constrain the sample set using the command line option"
                " '--samples'";
    }
    
    std::size_t num_samples_;
    
public:
    BadTrioSampleSet(std::size_t num_samples) : num_samples_ {num_samples} {}
};

class BadTrio : public UserError
{
    std::string do_where() const override
    {
        return "make_trio";
    }
    
    std::string do_why() const override
    {
        return "The given maternal and paternal samples are the same";
    }
    
    std::string do_help() const override
    {
        return "Ensure the sample names given in the command line options"
               " '--maternal-sample' and '--paternal-sample' differ and"
                " refer to valid samples";
    }
};

class BadTrioSamples : public UserError
{
    std::string do_where() const override
    {
        return "make_trio";
    }
    
    std::string do_why() const override
    {
        std::ostringstream ss {};
        if (mother_ && father_) {
            ss << "Neither of the parent sample names given command line options"
                  " '--maternal-sample' (" << *mother_ << ") and '--paternal-sample' ("
               << *father_ << ") appear in the read sample set";
        } else if (mother_) {
            ss << "The maternal sample name given in the command line option"
                    " '--maternal-sample' (" << *mother_ << ") does not appear in the"
                    " read sample set";
        } else {
            assert(father_);
            ss << "The paternal sample name given in the command line option"
            " '--paternal-sample' (" << *father_  << ") does not appear in the"
            " read sample set";
        }
        return ss.str();
    }
    
    std::string do_help() const override
    {
        return "Ensure the sample names given in the command line options"
        " '--maternal-sample' and '--paternal-sample' refer to valid samples";
    }
    
    boost::optional<SampleName> mother_, father_;
    
public:
    BadTrioSamples(boost::optional<SampleName> mother, boost::optional<SampleName> father)
    : mother_ {std::move(mother)}
    , father_ {std::move(father)}
    {}
};

auto get_caller_type(const OptionMap& options, const std::vector<SampleName>& samples,
	 				 const boost::optional<Pedigree>& pedigree)
{
    // TODO: could think about getting rid of the 'caller' option and just
    // deduce the caller type directly from the options.
    // Will need to report an error if conflicting caller options are given anyway.
    auto result = options.at("caller").as<std::string>();
    if (result == "population" && samples.size() == 1) {
        result = "individual";
    }
    if (options.count("maternal-sample") == 1 || options.count("paternal-sample") == 1
		|| (pedigree && is_trio(samples, *pedigree))) {
        result = "trio";
    }
    if (options.count("normal-sample") == 1) {
        result = "cancer";
    }
    return result;
}

auto get_child_from_trio(std::vector<SampleName> trio, const Pedigree& pedigree)
{
	if (is_parent_of(trio[0], trio[1], pedigree)) return trio[1];
	return is_parent_of(trio[1], trio[0], pedigree) ? trio[0] : trio[2];
}

Trio make_trio(std::vector<SampleName> samples, const Pedigree& pedigree)
{
	return *make_trio(get_child_from_trio(samples, pedigree), pedigree);
}

Trio make_trio(std::vector<SampleName> samples, const OptionMap& options,
			   const boost::optional<Pedigree>& pedigree)
{
	if (pedigree && is_trio(samples, *pedigree)) {
		return make_trio(samples, *pedigree);
	}
    if (samples.size() != 3) {
        throw BadTrioSampleSet {samples.size()};
    }
    auto mother = options.at("maternal-sample").as<SampleName>();
    auto father = options.at("paternal-sample").as<SampleName>();
    if (mother == father) {
        throw BadTrio {};
    }
    std::array<SampleName, 2> parents {mother, father};
    std::vector<SampleName> children {};
    std::sort(std::begin(samples), std::end(samples));
    std::sort(std::begin(parents), std::end(parents));
    assert(std::unique(std::begin(samples), std::end(samples)) == std::end(samples));
    std::set_difference(std::cbegin(samples), std::cend(samples),
                        std::cbegin(parents), std::cend(parents),
                        std::back_inserter(children));
    if (children.size() != 1) {
        boost::optional<SampleName> bad_mother, bad_father;
        if (!std::binary_search(std::cbegin(samples), std::cend(samples), mother)) {
            bad_mother = std::move(mother);
        }
        if (!std::binary_search(std::cbegin(samples), std::cend(samples), father)) {
            bad_father = std::move(father);
        }
        throw BadTrioSamples {std::move(bad_mother), std::move(bad_father)};
    }
    return Trio {
        Trio::Mother {std::move(mother)},
        Trio::Father {std::move(father)},
        Trio::Child  {std::move(children.front())}
    };
}

class UnimplementedCaller : public ProgramError
{
    std::string do_where() const override
    {
        return "get_caller_type";
    }
    
    std::string do_why() const override
    {
        return "The " + caller_ + " caller is not yet implemented. Sorry!";
    }
    
    std::string do_help() const override
    {
        return "please wait for updates";
    }
    
    std::string caller_;

public:
    UnimplementedCaller(std::string caller) : caller_ {caller} {}
};

bool allow_flank_scoring(const OptionMap& options)
{
    return options.at("inactive-flank-scoring").as<bool>() && !is_very_fast_mode(options);
}

CallerFactory make_caller_factory(const ReferenceGenome& reference, ReadPipe& read_pipe,
                                  const InputRegionMap& regions, const OptionMap& options)
{
    CallerBuilder vc_builder {
        reference,
        read_pipe,
        make_variant_generator_builder(options),
        make_haplotype_generator_builder(options)
    };
    
	const auto pedigree = get_pedigree(options);
    const auto caller = get_caller_type(options, read_pipe.samples(), pedigree);
    vc_builder.set_caller(caller);
    
    if (options.count("refcall") == 1) {
        const auto refcall_type = options.at("refcall").as<RefCallType>();
        if (refcall_type == RefCallType::positional) {
            vc_builder.set_refcall_type(CallerBuilder::RefCallType::positional);
        } else {
            vc_builder.set_refcall_type(CallerBuilder::RefCallType::blocked);
        }
        auto min_refcall_posterior = options.at("min-refcall-posterior").as<Phred<double>>();
        vc_builder.set_min_refcall_posterior(min_refcall_posterior);
    } else {
        vc_builder.set_refcall_type(CallerBuilder::RefCallType::none);
    }
    
    auto min_variant_posterior = options.at("min-variant-posterior").as<Phred<double>>();
    
    if (options.count("regenotype") == 1) {
        if (caller == "cancer") {
            vc_builder.set_min_variant_posterior(min_variant_posterior);
        } else {
            vc_builder.set_min_variant_posterior(Phred<double> {1});
        }
    } else {
        vc_builder.set_min_variant_posterior(min_variant_posterior);
    }
    vc_builder.set_ploidies(get_ploidy_map(options));
    vc_builder.set_max_haplotypes(get_max_haplotypes(options));
    vc_builder.set_haplotype_extension_threshold(options.at("haplotype-extension-threshold").as<Phred<double>>());
    auto min_phase_score = options.at("min-phase-score").as<Phred<double>>();
    vc_builder.set_min_phase_score(min_phase_score);
    if (!options.at("use-uniform-genotype-priors").as<bool>()) {
        vc_builder.set_snp_heterozygosity(options.at("snp-heterozygosity").as<float>());
        vc_builder.set_indel_heterozygosity(options.at("indel-heterozygosity").as<float>());
    }
    
    if (caller == "cancer") {
        if (options.count("normal-sample") == 1) {
            vc_builder.set_normal_sample(options.at("normal-sample").as<std::string>());
        }
        vc_builder.set_somatic_mutation_rate(options.at("somatic-mutation-rate").as<float>());
        vc_builder.set_min_expected_somatic_frequency(options.at("min-expected-somatic-frequency").as<float>());
        vc_builder.set_credible_mass(options.at("credible-mass").as<float>());
        vc_builder.set_min_credible_somatic_frequency(options.at("min-credible-somatic-frequency").as<float>());
        auto min_somatic_posterior = options.at("min-somatic-posterior").as<Phred<double>>();
        vc_builder.set_min_somatic_posterior(min_somatic_posterior);
    } else if (caller == "trio") {
        vc_builder.set_trio(make_trio(read_pipe.samples(), options, pedigree));
        vc_builder.set_snv_denovo_mutation_rate(options.at("snv-denovo-mutation-rate").as<float>());
        vc_builder.set_indel_denovo_mutation_rate(options.at("indel-denovo-mutation-rate").as<float>());
        vc_builder.set_min_denovo_posterior(options.at("min-denovo-posterior").as<Phred<double>>());
    }
    
    if (options.count("model-filtering") == 1) {
        vc_builder.set_model_filtering(options.at("model-filtering").as<bool>());
    } else {
        vc_builder.set_model_filtering(caller == "cancer");
    }
    
    if (call_sites_only(options)) {
        vc_builder.set_sites_only();
    }
    vc_builder.set_flank_scoring(allow_flank_scoring(options));
    vc_builder.set_model_mapping_quality(options.at("model-mapping-quality").as<bool>());
    vc_builder.set_max_joint_genotypes(as_unsigned("max-joint-genotypes", options));
    vc_builder.set_explain_read_directions(options.at("explain-read-directions").as<bool>());
    
    if (options.count("sequence-error-model") == 1) {
        vc_builder.set_sequencer(options.at("sequence-error-model").as<std::string>());
    }
    
    return CallerFactory {std::move(vc_builder)};
}

boost::optional<fs::path> get_output_path(const OptionMap& options)
{
    if (options.count("output") == 1) {
        return resolve_path(options.at("output").as<fs::path>(), options);
    }
    return boost::none;
}

VcfWriter make_output_vcf_writer(const OptionMap& options)
{
    auto output = get_output_path(options);
    return output ? VcfWriter {std::move(*output)} : VcfWriter {};
}

boost::optional<fs::path> create_temp_file_directory(const OptionMap& options)
{
    const auto working_directory = get_working_directory(options);
    auto result = working_directory;
    const fs::path temp_dir_base_name {"octopus-temp"};
    result /= temp_dir_base_name;
    constexpr unsigned temp_dir_name_count_limit {10000};
    unsigned temp_dir_counter {2};
    
    logging::WarningLogger log {};
    
    while (fs::exists(result) && temp_dir_counter <= temp_dir_name_count_limit) {
        if (fs::is_empty(result)) {
            stream(log) << "Found empty temporary directory " << result
            << ", it may need to be deleted manually";
        }
        result = working_directory;
        result /= temp_dir_base_name.string() + "-" + std::to_string(temp_dir_counter);
        ++temp_dir_counter;
    }
    
    if (temp_dir_counter > temp_dir_name_count_limit) {
        log << "There are many temporary directories in working directory indicating an error"
        " - new directory request blocked";
        return boost::none;
    }
    
    if (!fs::create_directory(result)) {
        stream(log) << "Failed to create temporary directory " << result << " - check permissions";
        return boost::none;
    }
    
    return result;
}

bool legacy_vcf_requested(const OptionMap& options)
{
    return options.at("legacy").as<bool>();
}
} // namespace options
} // namespace octopus

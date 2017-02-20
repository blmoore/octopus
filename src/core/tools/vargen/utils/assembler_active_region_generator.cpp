// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "assembler_active_region_generator.hpp"

#include <iterator>
#include <algorithm>
#include <cmath>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

#include "basics/genomic_region.hpp"
#include "basics/cigar_string.hpp"
#include "basics/aligned_read.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/maths.hpp"
#include "utils/append.hpp"

#include <iostream>

namespace octopus { namespace coretools {

AssemblerActiveRegionGenerator::AssemblerActiveRegionGenerator(const ReferenceGenome& reference)
: reference_ {reference}
, coverage_tracker_ {}
, interesting_read_coverages_ {}
{}

void AssemblerActiveRegionGenerator::add(const AlignedRead& read)
{
    coverage_tracker_.add(read);
    if (is_interesting(read)) {
        interesting_read_coverages_.add(read);
    }
}

template <typename Container>
void log_each(Container& values)
{
    for (auto& value : values) value = std::log(value);
}

template <typename Container>
auto to_logs(const Container& values)
{
    auto result = values;
    log_each(result);
    return result;
}

struct HiddenStateParameters
{
    double mean, stdev;
};

template <typename RealType>
RealType log_normal_pdf(const RealType x, const RealType mean, const RealType stdev)
{
    constexpr RealType c {0.9189385332046727417803297364056176398613974736377834};
    return std::log(stdev) - c - std::pow(x - mean, 2) / (2 * std::pow(stdev, 2));
}

double emmision_probability(const unsigned x, const HiddenStateParameters params)
{
    return log_normal_pdf(static_cast<double>(x), params.mean, params.stdev);
}

auto compute_base_deletion_probabilities(const std::vector<unsigned>& coverages,
                                         const double mean_coverage,
                                         const double stdev_coverage)
{
    constexpr std::size_t num_states {3};
    constexpr std::array<double, num_states> begin_probabilities {1 - 2e-10, 1e-10, 1e-10};
    constexpr std::array<std::array<double, num_states + 1>, num_states> transition_probabilities {{
                                                                                                   {1.0 - (1e-2 + 1e-3 + 1e-20), 1e-2, 1e-3, 1e-20}, // normal
                                                                                                   {1.0 - (1e-7 + 1e-15 + 1e-20), 1e-7, 1e-15, 1e-20}, // insertion
                                                                                                   {1e-5, 1e-15, 1.0 - (1e-5 + 1e-15 + 1e-20), 1e-20}, // deletion
                                                                                                   }};
    constexpr std::array<double, num_states> end_probabilities {1 - 2e-10, 1e-10, 1e-10};
    std::array<HiddenStateParameters, num_states> state_params {{{mean_coverage, stdev_coverage},
                                                                {mean_coverage + 3 * stdev_coverage, stdev_coverage},
                                                                {std::max(mean_coverage - 3 * stdev_coverage, 0.0), stdev_coverage / 2}}};
    const auto num_observations = coverages.size();
    // convert to log space
    const auto a_0 = to_logs(begin_probabilities);
    const auto a_e = to_logs(end_probabilities);
    auto a = transition_probabilities;
    for (auto& c : a) log_each(c);
    using State = std::array<double, num_states>;
    using maths::log_sum_exp;
    // Forward part
    std::vector<State> fwd(num_observations);
    fwd[0] = a_0;
    for (std::size_t s {0}; s < num_states; ++s) {
        fwd[0][s] += emmision_probability(coverages[0], state_params[s]);
    }
    for (std::size_t i {1}; i < num_observations; ++i) {
        for (std::size_t s {0}; s < num_states; ++s) {
            State tmp {};
            for (std::size_t t {0}; t < num_states; ++t) {
                tmp[t] = fwd[i - 1][t] + a[t][s];
            }
            fwd[i][s] = log_sum_exp(tmp) + emmision_probability(coverages[i], state_params[s]);
        }
    }
    auto end = fwd.back();
    std::transform(std::cbegin(end), std::cend(end), std::cbegin(a_e), std::begin(end),
                   [] (auto p, auto l) { return p + l; });
    const auto p_fwd = log_sum_exp(end);
    // Backward part
    std::vector<State> bkw(num_observations);
    bkw.back() = a_e;
    for (int i = num_observations - 2; i >= 0; --i) {
        for (std::size_t s {0}; s < num_states; ++s) {
            State tmp {};
            for (std::size_t t {0}; t < num_states; ++t) {
                tmp[t] = bkw[i + 1][t] + emmision_probability(coverages[i + 1], state_params[t]) + a[s][t];
            }
            bkw[i][s] = log_sum_exp(tmp);
        }
    }
    // Merge
    std::vector<State> posteriors(num_observations);
    std::transform(std::cbegin(fwd), std::cend(fwd), std::cbegin(bkw), std::begin(posteriors),
                   [=] (const auto f, const auto b) {
                       State result {};
                       for (std::size_t s {0}; s < num_states; ++s) {
                           result[s] = std::exp(f[s] + b[s] - p_fwd);
                       }
                       return result;
                   });
    std::vector<double> result(num_observations);
    std::transform(std::cbegin(posteriors), std::cend(posteriors), std::begin(result),
                   [] (const auto& s) { return s[2]; });
    return result;
}

auto get_regions(const std::vector<bool>& good_bases, const GenomicRegion& region)
{
    std::vector<GenomicRegion> result {};
    auto itr = std::find(std::cbegin(good_bases), std::cend(good_bases), true);
    for (; itr != std::cend(good_bases);) {
        const auto itr2 = std::find(itr, std::cend(good_bases), false);
        const auto begin = region.begin() + std::distance(std::cbegin(good_bases), itr);
        const auto end   = begin + std::distance(itr, itr2);
        result.emplace_back(region.contig_name(), begin, end);
        itr = std::find(itr2, std::cend(good_bases), true);
    }
    return result;
}

auto get_deletion_hotspots(const GenomicRegion& region, const CoverageTracker& tracker)
{
    const auto coverages = tracker.coverage(region);
    const auto mean_coverage = tracker.mean_coverage(region);
    const auto stdev_coverage = tracker.stdev_coverage(region);
    const auto deletion_base_probs = compute_base_deletion_probabilities(coverages, mean_coverage, stdev_coverage);
    std::vector<bool> deletion_bases(deletion_base_probs.size());
    std::transform(std::cbegin(deletion_base_probs), std::cend(deletion_base_probs), std::begin(deletion_bases),
                   [] (const auto p) {  return p > 0.5; });
    return get_regions(deletion_bases, region);
}

auto get_interesting_hotspots(const GenomicRegion& region, const CoverageTracker& interesting_read_tracker,
                              const CoverageTracker& tracker)
{
    const auto coverages = tracker.coverage(region);
    const auto interesting_coverages = interesting_read_tracker.coverage(region);
    std::vector<bool> interesting_bases(interesting_coverages.size());
    std::transform(std::cbegin(coverages), std::cend(coverages), std::cbegin(interesting_coverages),
                   std::begin(interesting_bases),
                   [] (const auto coverage, const auto interesting_coverage) {
                       return interesting_coverage >= coverage / 3;
                   });
    return get_regions(interesting_bases, region);
}

void merge(std::vector<GenomicRegion>&& src, std::vector<GenomicRegion>& dst)
{
    const auto itr = utils::append(std::move(src), dst);
    std::inplace_merge(std::begin(dst), itr, std::end(dst));
}

std::vector<GenomicRegion> AssemblerActiveRegionGenerator::generate(const GenomicRegion& region) const
{
    auto interesting_regions = get_interesting_hotspots(region, interesting_read_coverages_, coverage_tracker_);
    auto deletion_regions = get_deletion_hotspots(region, coverage_tracker_);
    merge(std::move(deletion_regions), interesting_regions);
    return join(extract_covered_regions(interesting_regions), 30);
}

// private methods

namespace {

using NucleotideSequenceIterator = AlignedRead::NucleotideSequence::const_iterator;
using BaseQualityVectorIterator = AlignedRead::BaseQualityVector::const_iterator;

bool has_snv_in_match_range(const NucleotideSequenceIterator first_ref, const NucleotideSequenceIterator last_ref,
                            const NucleotideSequenceIterator first_base, const BaseQualityVectorIterator first_quality,
                            const AlignedRead::BaseQuality trigger = 20)
{
    using boost::make_zip_iterator;
    using Tuple = boost::tuple<char, char, AlignedRead::BaseQuality>;
    const auto num_bases = std::distance(first_ref, last_ref);
    const auto last_base = std::next(first_base, num_bases);
    const auto last_quality = std::next(first_quality, num_bases);
    return std::any_of(make_zip_iterator(boost::make_tuple(first_ref, first_base, first_quality)),
                       make_zip_iterator(boost::make_tuple(last_ref, last_base, last_quality)),
                       [trigger] (const Tuple& t) {
                           const char ref_base  {t.get<0>()}, read_base {t.get<1>()};
                           return ref_base != read_base && ref_base != 'N' && read_base != 'N' && t.get<2>() >= trigger;
                       });
}

bool is_good_clip(const NucleotideSequenceIterator first_base, const NucleotideSequenceIterator last_base,
                  const BaseQualityVectorIterator first_quality, const AlignedRead::BaseQuality good_base_trigger,
                  const std::size_t min_good_bases)
{
    using boost::make_zip_iterator;
    using Tuple = boost::tuple<char, AlignedRead::BaseQuality>;
    const auto last_quality = std::next(first_quality, std::distance(first_base, last_base));
    return static_cast<std::size_t>(std::count_if(make_zip_iterator(boost::make_tuple(first_base, first_quality)),
                                                  make_zip_iterator(boost::make_tuple(last_base, last_quality)),
                                                  [=] (const Tuple& t) {
                                                      return t.get<0>() != 'N' && t.get<1>() >= good_base_trigger;
                                                  })) >= min_good_bases;
}
    
} // namespace

bool AssemblerActiveRegionGenerator::is_interesting(const AlignedRead& read) const
{
    using std::cbegin; using std::next; using std::move;
    using Flag = CigarOperation::Flag;
    const auto& read_sequence = read.sequence();
    auto sequence_itr = cbegin(read_sequence);
    auto base_quality_itr = cbegin(read.qualities());
    auto ref_index = mapped_begin(read);
    std::size_t read_index {0};
    constexpr AlignedRead::BaseQuality trigger_quality {20};
    constexpr CigarOperation::Size trigger_clip_size {5};
    for (const auto& cigar_operation : read.cigar()) {
        const auto op_size = cigar_operation.size();
        switch (cigar_operation.flag()) {
            case Flag::alignmentMatch:
            {
//                const GenomicRegion region{contig_name(read), ref_index, ref_index + op_size};
//                const auto ref_segment = reference_.get().fetch_sequence(region);
//                if (has_snv_in_match_range(std::cbegin(ref_segment), std::cend(ref_segment),
//                                           next(sequence_itr, read_index),
//                                           next(base_quality_itr, read_index),
//                                           trigger_quality)) {
//                    return true;
//                }
                read_index += op_size;
                ref_index += op_size;
                break;
            }
            case Flag::sequenceMatch:
                read_index += op_size;
                ref_index += op_size;
                break;
            case Flag::substitution:
            {
//                GenomicRegion {contig_name(read), ref_index, ref_index + op_size};
//                if (std::any_of(next(base_quality_itr, read_index), next(base_quality_itr, read_index + op_size),
//                                [trigger_quality] (const auto quality) { return quality >= trigger_quality; })) {
//                    return true;
//                }
                read_index += op_size;
                ref_index += op_size;
                break;
            }
            case Flag::insertion: return true;
            case Flag::deletion: return true;
            case Flag::softClipped:
            {
                if (is_good_clip(next(sequence_itr, read_index), next(sequence_itr, read_index + op_size),
                                 next(base_quality_itr, read_index), trigger_quality, trigger_clip_size)) {
                    return true;
                }
                read_index += op_size;
                ref_index += op_size;
                break;
            }
            case Flag::hardClipped: break;
            case Flag::padding:
                ref_index += op_size;
                break;
            case Flag::skipped:
                ref_index += op_size;
                break;
        }
    }
    return false;
}
    
} // namespace coretools
} // namespace octopus

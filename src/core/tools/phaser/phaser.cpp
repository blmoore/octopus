// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "phaser.hpp"

#include <deque>
#include <algorithm>
#include <numeric>
#include <cstddef>
#include <cmath>
#include <utility>
#include <iostream>

#include <boost/functional/hash.hpp>

#include "utils/mappable_algorithms.hpp"
#include "utils/maths.hpp"

#include "timers.hpp"

namespace octopus {

Phaser::Phaser(Phred<double> min_phase_score) : min_phase_score_ {min_phase_score} {}

namespace {

using GenotypeReference    = std::reference_wrapper<const Genotype<Haplotype>>;
using PhaseComplementSet   = std::deque<GenotypeReference>;
using PhaseComplementSets  = std::vector<PhaseComplementSet>;
using PartitionIterator    = std::vector<GenomicRegion>::const_iterator;
using GenotypeSpliceVector = std::vector<Genotype<Haplotype>>;

struct GenotypePartitionSplices
{
    GenotypeReference genotype;
    GenotypeSpliceVector splices;
};

using SpliceTable = std::vector<GenotypePartitionSplices>;

template <typename Container>
SpliceTable make_splice_table(const Container& genotypes, PartitionIterator first_parition, PartitionIterator last_partition)
{
    const auto num_partitions = std::distance(first_parition, last_partition);
    SpliceTable result {};
    result.reserve(genotypes.size());
    for (const auto& genotype : genotypes) {
        result.push_back({std::cref(genotype), GenotypeSpliceVector {}});
        result.back().splices.reserve(num_partitions);
    }
    std::for_each(first_parition, last_partition, [&] (const auto& region) {
        auto splices = splice_each<Haplotype>(genotypes, region);
        for (std::size_t i {0}; i < genotypes.size(); ++i) {
            result[i].splices.push_back(std::move(splices[i]));
        }
    });
    return result;
}

struct GenotypeSpliceVectorHash
{
    auto operator()(const GenotypeSpliceVector& splices) const noexcept
    {
        std::size_t result {};
        for (const auto& genotype : splices) {
            using boost::hash_combine;
            hash_combine(result, std::hash<Genotype<Haplotype>>{}(genotype));
        }
        return result;
    }
};

using PhaseSetMap = std::unordered_map<GenotypeSpliceVector, PhaseComplementSet, GenotypeSpliceVectorHash>;

PhaseComplementSets make_phase_sets(const SpliceTable& splice_table)
{
    PhaseSetMap phase_sets {};
    phase_sets.reserve(splice_table.size());
    for (const auto& row : splice_table) {
        phase_sets[row.splices].push_back(row.genotype);
    }
    PhaseComplementSets result {};
    result.reserve(phase_sets.size());
    for (auto&& p : phase_sets) {
        result.push_back(std::move(p.second));
    }
    return result;
}

template <typename Container>
PhaseComplementSets
generate_phase_complement_sets(const Container& genotypes, PartitionIterator first_parition, PartitionIterator last_partition)
{
    const auto splice_table = make_splice_table(genotypes, first_parition, last_partition);
    return make_phase_sets(splice_table);
}

template <typename Map>
double marginalise(const PhaseComplementSet& phase_set, const Map& genotype_posteriors)
{
    return std::accumulate(std::cbegin(phase_set), std::cend(phase_set), 0.0,
                           [&] (const auto curr, const auto& genotype) {
                               return curr + genotype_posteriors.at(genotype);
                           });
}

template <typename Map>
double calculate_entropy(const PhaseComplementSet& phase_set, const Map& genotype_posteriors)
{
    const auto norm = marginalise(phase_set, genotype_posteriors);
    return std::max(0.0, -std::accumulate(std::cbegin(phase_set), std::cend(phase_set), 0.0,
                                          [&genotype_posteriors, norm] (const auto curr, const auto& genotype) {
                                              const auto p = genotype_posteriors.at(genotype) / norm;
                                              return curr + p * std::log2(p);
                                          }));
}

double maximum_entropy(const size_t num_elements)
{
    return std::log2(num_elements);
}

template <typename Map>
double calculate_relative_entropy(const PhaseComplementSet& phase_set, const Map& genotype_posteriors)
{
    if (phase_set.size() < 2) return 1.0;
    return 1.0 - calculate_entropy(phase_set, genotype_posteriors) / maximum_entropy(phase_set.size());
}

template <typename Map>
auto calculate_phase_score(const PhaseComplementSet& phase_set, const Map& genotype_posteriors)
{
    return marginalise(phase_set, genotype_posteriors) * calculate_relative_entropy(phase_set, genotype_posteriors);
}

template <typename Map>
Phred<double> calculate_phase_score(const PhaseComplementSets& phase_sets, const Map& genotype_posteriors)
{
    return Phred<double> { Phred<double>::Probability {
        std::max(0.0, 1.0 - std::accumulate(std::cbegin(phase_sets), std::cend(phase_sets), 0.0,
                                           [&] (const auto curr, const auto& phase_set) {
                                               return curr + calculate_phase_score(phase_set, genotype_posteriors);
                                           }))
    }};
}

auto min_phase_score(const double p, const std::size_t num_genotypes)
{
    if (maths::almost_one(p) || num_genotypes == 1) {
        static const Phred<double> max_possible_score {Phred<double>::Probability {0.0}};
        return max_possible_score;
    } else {
        const auto s = p * std::log2(p) + (1.0 - p) * std::log2((1.0 - p) / (num_genotypes - 1));
        const auto y = 1.0 + s / maximum_entropy(2);
        return Phred<double> {Phred<double>::Probability {std::max(y, 0.0)}};
    }
}

auto min_phase_score(const Genotype<Haplotype>& called_genotype, const Phaser::SampleGenotypePosteriorMap& genotype_posteriors)
{
    return min_phase_score(genotype_posteriors[called_genotype], genotype_posteriors.size());
}

std::vector<GenotypeReference> extract_genotypes(const Phaser::GenotypePosteriorMap& genotype_posteriors)
{
    return extract_key_refs(genotype_posteriors);
}

using GenotypeSplicePosteriorMap = std::unordered_map<Genotype<Haplotype>, double>;

template <typename Container>
auto splice_and_marginalise(const Container& genotypes,
                            const Phaser::SampleGenotypePosteriorMap& genotype_posteriors,
                            const GenomicRegion& region)
{
    auto splices = splice_each<Haplotype>(genotypes, region);
    GenotypeSplicePosteriorMap splice_posteriors {genotypes.size()};
    auto splice_itr = std::cbegin(splices);
    for (const auto& genotype : genotypes) {
        splice_posteriors[*splice_itr++] += genotype_posteriors[genotype];
    }
    std::sort(std::begin(splices), std::end(splices), GenotypeLess {});
    splices.erase(std::unique(std::begin(splices), std::end(splices)), std::end(splices));
    return std::make_pair(std::move(splices), std::move(splice_posteriors));
}

} // namespace

boost::optional<Phaser::PhaseSet>
Phaser::try_phase(const std::vector<Haplotype>& haplotypes,
                  const GenotypePosteriorMap& genotype_posteriors,
                  const std::vector<GenomicRegion>& regions) const
{
    // TODO
    return boost::none;
}

Phaser::PhaseSet::SamplePhaseRegions
force_phase_sample(const GenomicRegion& region,
                   const std::vector<GenomicRegion>& partitions,
                   const std::vector<GenotypeReference>& genotypes,
                   const Phaser::SampleGenotypePosteriorMap& genotype_posteriors,
                   const Phred<double> min_phase_score)
{
    auto first_partition = std::cbegin(partitions);
    auto last_partition  = std::cend(partitions);
    auto phase_set = generate_phase_complement_sets(genotypes, first_partition, last_partition);
    auto phase_score = calculate_phase_score(phase_set, genotype_posteriors);
    if (phase_score >= min_phase_score) {
        return {Phaser::PhaseSet::PhaseRegion {region, phase_score}};
    }
    Phaser::PhaseSet::SamplePhaseRegions result {};
    --last_partition;
    std::vector<Genotype<Haplotype>> splices;
    GenotypeSplicePosteriorMap splice_posteriors;
    while (first_partition != std::cend(partitions)) {
        auto curr_region = encompassing_region(first_partition, last_partition);
        std::tie(splices, splice_posteriors) = splice_and_marginalise(genotypes, genotype_posteriors, curr_region);
        phase_set = generate_phase_complement_sets(splices, first_partition, last_partition);
        phase_score = calculate_phase_score(phase_set, splice_posteriors);
        if (phase_score >= min_phase_score || std::distance(first_partition, last_partition) == 1) {
            result.emplace_back(encompassing_region(first_partition, last_partition), phase_score);
            first_partition = last_partition;
            last_partition  = std::cend(partitions);
        } else {
            --last_partition;
        }
    }
    return result;
}

Phaser::PhaseSet
Phaser::force_phase(const std::vector<Haplotype>& haplotypes,
                    const GenotypePosteriorMap& genotype_posteriors,
                    const std::vector<GenomicRegion>& regions,
                    boost::optional<GenotypeCallMap> genotype_calls) const
{
    assert(!haplotypes.empty());
    assert(!genotype_posteriors.empty1());
    assert(!genotype_posteriors.empty2());
    assert(std::is_sorted(std::cbegin(regions), std::cend(regions)));
    const auto& haplotype_region = haplotypes.front().mapped_region();
    const auto genotypes = extract_genotypes(genotype_posteriors);
    const auto partitions = extract_covered_regions(regions);
    PhaseSet result {haplotype_region};
    result.phase_regions.reserve(genotype_posteriors.size1());
    if (genotypes.front().get().ploidy() == 1 || partitions.size() == 1) {
        for (const auto& p : genotype_posteriors) {
            if (max_phase_score_) {
                result.phase_regions[p.first].emplace_back(haplotype_region, *max_phase_score_);
            } else {
                static const Phred<double> max_possible_score {Phred<double>::Probability {0.0}};
                result.phase_regions[p.first].emplace_back(haplotype_region, max_possible_score);
            }
        }
        return result;
    }
    for (const auto& p : genotype_posteriors) {
        if (genotype_calls && max_phase_score_ && min_phase_score(genotype_calls->at(p.first), p.second) >= *max_phase_score_) {
            result.phase_regions[p.first].emplace_back(haplotype_region, *max_phase_score_);
        } else {
            auto phases = force_phase_sample(haplotype_region, partitions, genotypes, p.second, min_phase_score_);
            if (max_phase_score_) {
                for (auto& phase : phases) phase.score = std::min(phase.score, *max_phase_score_);
            }
            result.phase_regions.emplace(p.first, std::move(phases));
        }
    }
    return result;
}

// non-member methods

bool is_split_phasing(const Phaser::PhaseSet& phase)
{
    return std::any_of(std::cbegin(phase.phase_regions), std::cend(phase.phase_regions),
                       [] (const auto& p) { return p.second.size() > 1; });
}
    
namespace debug {

void print_phase_sets(const Phaser::PhaseSet& phasings)
{
    print_phase_sets(std::cout, phasings);
}

} // namespace debug

} // namespace Ocotpus

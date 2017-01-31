// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_direction_model.hpp"

#include <cmath>
#include <cstddef>

#include "utils/maths.hpp"

namespace octopus {

double log_binomial(const unsigned n, const unsigned k)
{
    using maths::log_factorial;
    return log_factorial<double>(n) - (log_factorial<double>(n - k) + log_factorial<double>(k));
}

double calculate_read_direction_probability(const Genotype<Haplotype>& genotype,
                                            const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                            const std::vector<AlignedRead::Direction>& read_directions)
{
    if (genotype.is_homozygous()) {
        //const auto num_forward = std::accumulate(std::cbegin(read_directions), std::cend(read_directions), 0);
        //return log_binomial(static_cast<unsigned>(read_directions.size()), static_cast<unsigned>(num_forward));
        return 0;
    }
    double hap1_forward {0}, hap1_reverse {0}, hap2_forward {}, hap2_reverse {0};
    const auto& hap1_likelihoods = haplotype_likelihoods[genotype[0]];
    const auto& hap2_likelihoods = haplotype_likelihoods[genotype[1]];
    for (std::size_t i {0}; i < read_directions.size(); ++i) {
        if (hap1_likelihoods[i] > hap2_likelihoods[i]) {
            if (read_directions[i] == AlignedRead::Direction::forward) {
                ++hap1_forward;
            } else {
                ++hap1_reverse;
            }
        } else if (maths::almost_equal(hap1_likelihoods[i], hap2_likelihoods[i])) {
//            if (read_directions[i]) {
//                hap1_forward += 0.5;
//                hap2_forward += 0.5;
//            } else {
//                hap1_reverse += 0.5;
//                hap2_reverse += 0.5;
//            }
        } else {
            if (read_directions[i] == AlignedRead::Direction::forward) {
                ++hap2_forward;
            } else {
                ++hap2_reverse;
            }
        }
    }
    const auto num_hap1 = static_cast<unsigned>(hap1_forward + hap1_reverse);
    const auto num_hap2 = static_cast<unsigned>(hap2_forward + hap2_reverse);
    const auto b1 = log_binomial(num_hap1, static_cast<unsigned>(hap1_forward));
    const auto b2 = log_binomial(num_hap2, static_cast<unsigned>(hap2_forward));
    return (b1 + b2) - (std::log(2) * (num_hap1 + num_hap2));
}

} // namespace octopus

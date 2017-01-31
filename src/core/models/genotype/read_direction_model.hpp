// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_direction_model_hpp
#define read_direction_model_hpp

#include <vector>

#include "basics/aligned_read.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/models/haplotype_likelihood_cache.hpp"

namespace octopus {

// p(read_directions | genotype, haplotype_likelihoods)
double calculate_read_direction_probability(const Genotype<Haplotype>& genotype,
                                            const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                            const std::vector<AlignedRead::Direction>& read_directions);

} // namespace octopus


#endif

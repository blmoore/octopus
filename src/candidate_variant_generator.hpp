//
//  candidate_variant_generator.hpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__candidate_variant_generator__
#define __Octopus__candidate_variant_generator__

#include <vector>
#include <memory>    // std::unique_ptr
#include <cstddef>   // size_t
#include <iterator>  // std::cbegin, std::cend

#include "i_candidate_variant_generator.hpp"

class GenomicRegion;
class AlignedRead;
class Variant;

class CandidateVariantGenerator : public ICandidateVariantGenerator
{
public:
    CandidateVariantGenerator()           = default;
    ~CandidateVariantGenerator() override = default;
    
    CandidateVariantGenerator(const CandidateVariantGenerator&)            = default;
    CandidateVariantGenerator& operator=(const CandidateVariantGenerator&) = default;
    CandidateVariantGenerator(CandidateVariantGenerator&&)                 = default;
    CandidateVariantGenerator& operator=(CandidateVariantGenerator&&)      = default;
    
    void register_generator(std::unique_ptr<ICandidateVariantGenerator> generator);
    void add_read(const AlignedRead& a_read) override;
    void add_reads(std::vector<AlignedRead>::const_iterator first, std::vector<AlignedRead>::const_iterator last) override;
    void add_reads(MappableSet<AlignedRead>::const_iterator first, MappableSet<AlignedRead>::const_iterator last) override;
    std::vector<Variant> get_candidates(const GenomicRegion& a_region) override;
    void reserve(size_t n) override;
    void clear() override;
    
private:
    std::vector<std::unique_ptr<ICandidateVariantGenerator>> generator_list_;
};

template <typename Map>
void add_reads(const Map& map, CandidateVariantGenerator& generator)
{
    for (const auto& m : map) {
        generator.add_reads(std::cbegin(m.second), std::cend(m.second));
    }
}

#endif /* defined(__Octopus__variant_candidate_generator__) */
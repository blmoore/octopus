// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "mapping_quality_zero_count.hpp"

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "utils/read_stats.hpp"
#include "../facets/overlapping_reads.hpp"

namespace octopus { namespace csr {

MappingQualityZeroCount::MappingQualityZeroCount(bool recalculate) : recalculate_ {recalculate} {}

std::unique_ptr<Measure> MappingQualityZeroCount::do_clone() const
{
    return std::make_unique<MappingQualityZeroCount>(*this);
}

Measure::ResultType MappingQualityZeroCount::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    if (recalculate_) {
        const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
        return count_mapq_zero(reads);
    } else {
        return static_cast<std::size_t>(std::stoull(call.info_value("MQ0").front()));
    }
}

std::string MappingQualityZeroCount::do_name() const
{
    return "MQ0";
}

std::vector<std::string> MappingQualityZeroCount::do_requirements() const
{
    if (recalculate_) {
        return {"OverlappingReads"};
    } else {
        return {};
    }
}

} // namespace csr
} // namespace octopus

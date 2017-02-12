// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "caller_wrapper.hpp"

#include <utility>

#include "core/callers/caller.hpp"
#include "logging/progress_meter.hpp"

namespace octopus {

CallerWrapper::CallerWrapper(std::unique_ptr<Caller> caller, ProgressMeter& progress_meter,
                             boost::optional<VariantGenerator> regenotype_variant_generator)
: caller_ {std::move(caller)}
, progress_meter_ {progress_meter}
, regenotype_variant_generator_ {std::move(regenotype_variant_generator)}
{}

std::vector<VcfRecord> CallerWrapper::call(const GenomicRegion& region) const
{
	if (regenotype_variant_generator_) {
        const auto regenotype_variants = regenotype_variant_generator_->generate(region);
        return caller_->regenotype(regenotype_variants, progress_meter_);
    } else {
        return caller_->call(region, progress_meter_);
    }
}

} // namespace octopus

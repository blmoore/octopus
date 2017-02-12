// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <memory>
#include <vector>

#include <boost/optional.hpp>

#include "core/tools/vargen/variant_generator.hpp"

namespace octopus {

class Caller;
class VcfRecord;
class GenomicRegion;
class ProgressMeter;

class CallerWrapper
{
public:
    CallerWrapper() = delete;
    
    CallerWrapper(std::unique_ptr<Caller> caller, ProgressMeter& progress_meter,
		 		  boost::optional<VariantGenerator> regenotype_variant_generator = boost::none);
	
    CallerWrapper(const CallerWrapper&)             = delete;
    CallerWrapper& operator=(const CallerWrapper&)  = delete;
    CallerWrapper(CallerWrapper&& other) 		    = default;
    CallerWrapper& operator=(CallerWrapper&& other) = default;
    
    ~CallerWrapper() = default;
	
	std::vector<VcfRecord> call(const GenomicRegion& region) const;
	
private:
	std::unique_ptr<Caller> caller_;
	ProgressMeter& progress_meter_;
	mutable boost::optional<VariantGenerator> regenotype_variant_generator_;
};

} // namespace octopus

// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variant_call_filter_hpp
#define variant_call_filter_hpp

#include <vector>
#include <string>
#include <cstddef>
#include <type_traits>
#include <functional>
#include <future>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/phred.hpp"
#include "basics/genomic_region.hpp"
#include "core/types/variant.hpp"
#include "io/variant/vcf_header.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_reader.hpp"
#include "utils/thread_pool.hpp"
#include "logging/logging.hpp"
#include "../facets/facet.hpp"
#include "../facets/facet_factory.hpp"
#include "../measures/measure.hpp"

namespace octopus {

class VcfWriter;
class VcfHeader;

namespace csr {

class VariantCallFilter
{
public:
    struct OutputOptions
    {
        bool emit_sites_only = false;
        bool clear_existing_filters = true;
        bool annotate_measures = false;
        bool clear_info = false;
    };
    
    struct ConcurrencyPolicy
    {
        boost::optional<unsigned> max_threads = boost::none;
    };
    
    VariantCallFilter() = delete;
    
    VariantCallFilter(FacetFactory facet_factory,
                      std::vector<MeasureWrapper> measures,
                      OutputOptions output_config,
                      ConcurrencyPolicy threading);
    
    VariantCallFilter(const VariantCallFilter&)            = delete;
    VariantCallFilter& operator=(const VariantCallFilter&) = delete;
    VariantCallFilter(VariantCallFilter&&)                 = default;
    VariantCallFilter& operator=(VariantCallFilter&&)      = default;
    
    virtual ~VariantCallFilter() = default;
    
    void filter(const VcfReader& source, VcfWriter& dest) const;
    
protected:
    using SampleList    = std::vector<SampleName>;
    using MeasureVector = std::vector<Measure::ResultType>;
    using VcfIterator   = VcfReader::RecordIterator;
    using CallBlock     = std::vector<VcfRecord>;
    using MeasureBlock  = std::vector<MeasureVector>;
    
    struct Classification
    {
        enum class Category { hard_filtered, soft_filtered, unfiltered } category;
        std::vector<std::string> reasons = {};
        boost::optional<Phred<double>> quality = boost::none;
    };
    
    mutable boost::optional<logging::DebugLogger> debug_log_;
    
    bool can_measure_single_call() const noexcept;
    bool can_measure_multiple_blocks() const noexcept;
    CallBlock read_next_block(VcfIterator& first, const VcfIterator& last, const SampleList& samples) const;
    std::vector<CallBlock> read_next_blocks(VcfIterator& first, const VcfIterator& last, const SampleList& samples) const;
    MeasureVector measure(const VcfRecord& call) const;
    MeasureBlock measure(const CallBlock& block) const;
    std::vector<MeasureBlock> measure(const std::vector<CallBlock>& blocks) const;
    void write(const VcfRecord& call, const Classification& classification, VcfWriter& dest) const;
    void annotate(VcfRecord::Builder& call, const MeasureVector& measures) const;
    
private:
    using FacetNameSet = std::vector<std::string>;
    
    FacetFactory facet_factory_;
    FacetNameSet facet_names_;
    std::vector<MeasureWrapper> measures_;
    OutputOptions output_config_;
    
    mutable ThreadPool workers_;
    
    virtual void annotate(VcfHeader::Builder& header) const = 0;
    virtual void filter(const VcfReader& source, VcfWriter& dest, const SampleList& samples) const = 0;
    
    VcfHeader make_header(const VcfReader& source) const;
    Measure::FacetMap compute_facets(const CallBlock& block) const;
    std::vector<Measure::FacetMap> compute_facets(const std::vector<CallBlock>& blocks) const;
    MeasureBlock measure(const CallBlock& block, const Measure::FacetMap& facets) const;
    MeasureVector measure(const VcfRecord& call, const Measure::FacetMap& facets) const;
    VcfRecord::Builder construct_template(const VcfRecord& call) const;
    void annotate(VcfRecord::Builder& call, Classification status) const;
    void pass(VcfRecord::Builder& call) const;
    void fail(VcfRecord::Builder& call, std::vector<std::string> reasons) const;
    bool is_multithreaded() const noexcept;
    unsigned max_concurrent_blocks() const noexcept;
};

} // namespace csr

using csr::VariantCallFilter;

} // namespace octopus

#endif

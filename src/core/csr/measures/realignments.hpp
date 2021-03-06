// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef realignments_hpp
#define realignments_hpp

#include <string>
#include <vector>

#include "measure.hpp"

namespace octopus {

class VcfRecord;

namespace csr {

class Realignments : public Measure
{
    std::unique_ptr<Measure> do_clone() const override;
    ResultType do_evaluate(const VcfRecord& call, const FacetMap& facets) const override;
    std::string do_name() const override;
    std::vector<std::string> do_requirements() const override;
    std::string do_serialise(const ResultType& value) const override;
};

} // namespace csr
} // namespace octopus

#endif
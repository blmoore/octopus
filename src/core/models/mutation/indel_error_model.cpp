// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "indel_error_model.hpp"

#include "core/types/haplotype.hpp"

namespace octopus {

void IndelErrorModel::set_penalities(const Haplotype& haplotype, PenaltyVector& gap_open, PenaltyVector& gap_extend) const
{
    return do_set_penalities(haplotype, gap_open, gap_extend);
}

} // namespace octopus

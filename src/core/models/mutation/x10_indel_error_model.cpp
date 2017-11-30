// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "x10_indel_error_model.hpp"

#include <algorithm>
#include <iterator>

#include "tandem/tandem.hpp"

#include "core/types/haplotype.hpp"

namespace octopus {

constexpr decltype(X10IndelErrorModel::homopolymerErrors_) X10IndelErrorModel::homopolymerErrors_;
constexpr decltype(X10IndelErrorModel::homopolymerErrors_) X10IndelErrorModel::diNucleotideTandemRepeatErrors_;
constexpr decltype(X10IndelErrorModel::homopolymerErrors_) X10IndelErrorModel::triNucleotideTandemRepeatErrors_;
constexpr decltype(X10IndelErrorModel::homopolymerErrors_) X10IndelErrorModel::polyNucleotideTandemRepeatErrors_;
constexpr decltype(X10IndelErrorModel::extendPenalties_) X10IndelErrorModel::extendPenalties_;

namespace {

auto extract_repeats(const Haplotype& haplotype)
{
    return tandem::extract_exact_tandem_repeats(haplotype.sequence(), 1, 3);
}

template <typename C, typename T>
static auto get_penalty(const C& penalties, const T length) noexcept
{
    return (length < penalties.size()) ? penalties[length] : penalties.back();
}

} // namespace

void X10IndelErrorModel::do_set_penalities(const Haplotype& haplotype, PenaltyVector& gap_open, PenaltyVector& gap_extend) const
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::next;
    const auto repeats = extract_repeats(haplotype);
    gap_open.assign(sequence_size(haplotype), homopolymerErrors_.front());
    for (const auto& repeat : repeats) {
        PenaltyType e;
        switch (repeat.period) {
            case 1:
            {
                e = get_penalty(homopolymerErrors_, repeat.length);
                break;
            }
            case 2:
            {
                static constexpr std::array<char, 2> AC {'A', 'C'};
                e = get_penalty(diNucleotideTandemRepeatErrors_, repeat.length / 2);
                const auto it = next(cbegin(haplotype.sequence()), repeat.pos);
                if (e > 10 && std::equal(cbegin(AC), cend(AC), it)) {
                    e -= 2;
                }
                break;
            }
            case 3:
            {
                static constexpr std::array<char, 3> GGC {'G', 'G', 'C'};
                static constexpr std::array<char, 3> GCC {'G', 'C', 'C'};
                e = get_penalty(triNucleotideTandemRepeatErrors_, repeat.length / 3);
                const auto it = next(cbegin(haplotype.sequence()), repeat.pos);
                if (e > 10 && std::equal(cbegin(GGC), cend(GGC), it)) {
                    e -= 2;
                } else if (e > 12 && std::equal(cbegin(GCC), cend(GCC), it)) {
                    e -= 1;
                }
                break;
            }
            default:
                e = get_penalty(polyNucleotideTandemRepeatErrors_, repeat.length / repeat.period);
        }
        std::fill_n(next(begin(gap_open), repeat.pos), repeat.length, e);
        auto local_gap_extend = get_penalty(extendPenalties_, repeat.length / repeat.period);
        std::fill_n(next(begin(gap_extend), repeat.pos), repeat.length, local_gap_extend);
    }
}
    
} // namespace octopus

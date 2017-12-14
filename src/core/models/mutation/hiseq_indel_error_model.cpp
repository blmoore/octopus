// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "hiseq_indel_error_model.hpp"

#include <algorithm>
#include <iterator>

#include "tandem/tandem.hpp"

#include "core/types/haplotype.hpp"

namespace octopus {

constexpr decltype(HiSeqIndelErrorModel::homopolymerErrors_) HiSeqIndelErrorModel::homopolymerErrors_;
constexpr decltype(HiSeqIndelErrorModel::homopolymerErrors_) HiSeqIndelErrorModel::diNucleotideTandemRepeatErrors_;
constexpr decltype(HiSeqIndelErrorModel::homopolymerErrors_) HiSeqIndelErrorModel::triNucleotideTandemRepeatErrors_;
constexpr decltype(HiSeqIndelErrorModel::homopolymerErrors_) HiSeqIndelErrorModel::polyNucleotideTandemRepeatErrors_;
constexpr decltype(HiSeqIndelErrorModel::extendPenalties_) HiSeqIndelErrorModel::extendPenalties_;

namespace {

auto extract_repeats(const Haplotype& haplotype)
{
    return tandem::extract_exact_tandem_repeats(haplotype.sequence(), 1, 3);
}

template <typename C, typename T>
static auto get_penalty(const C& penalties, const T length) noexcept
{
    return length < penalties.size() ? penalties[length] : penalties.back();
}

template <typename FordwardIt, typename Tp>
auto fill_if_less(FordwardIt first, FordwardIt last, const Tp& value)
{
    return std::transform(first, last, first, [&] (const auto& x) { return std::min(x, value); });
}

template <typename FordwardIt, typename Tp>
auto fill_n_if_less(FordwardIt first, std::size_t n, const Tp& value)
{
    return fill_if_less(first, std::next(first, n), value);
}

} // namespace

HiSeqIndelErrorModel::PenaltyType
HiSeqIndelErrorModel::do_evaluate(const Haplotype& haplotype, PenaltyVector& gap_open_penalities) const
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::next;
    const auto repeats = extract_repeats(haplotype);
    gap_open.assign(sequence_size(haplotype), homopolymerErrors_.front());
    gap_extend.assign(sequence_size(haplotype), extendPenalties_.front());
    for (const auto& repeat : repeats) {
        PenaltyType local_gap_open {};
        switch (repeat.period) {
            case 1:
            {
                local_gap_open = get_penalty(homopolymerErrors_, repeat.length);
                break;
            }
            case 2:
            {
                static constexpr std::array<char, 2> AC {'A', 'C'};
                local_gap_open = get_penalty(diNucleotideTandemRepeatErrors_, repeat.length / 2);
                if (local_gap_open > 10 && std::equal(cbegin(AC), cend(AC), next(cbegin(haplotype.sequence()), repeat.pos))) {
                    local_gap_open -= 2;
                }
                break;
            }
            case 3:
            {
                static constexpr std::array<char, 3> GGC {'G', 'G', 'C'};
                static constexpr std::array<char, 3> GCC {'G', 'C', 'C'};
                local_gap_open = get_penalty(triNucleotideTandemRepeatErrors_, repeat.length / 3);
                const auto itr = next(cbegin(haplotype.sequence()), repeat.pos);
                if (local_gap_open > 10 && std::equal(cbegin(GGC), cend(GGC), itr)) {
                    local_gap_open -= 2;
                } else if (local_gap_open > 12 && std::equal(cbegin(GCC), cend(GCC), itr)) {
                    local_gap_open -= 1;
                }
                break;
            }
            default:
                local_gap_open = get_penalty(polyNucleotideTandemRepeatErrors_, repeat.length / repeat.period);
        }
        default:
            local_gap_open = get_penalty(polyNucleotideTandemRepeatErrors_, repeat.length / repeat.period);
        }
        fill_n_if_less(next(begin(gap_open_penalities), repeat.pos), repeat.length, local_gap_open);
        if (repeat.length > max_repeat.length) {
            max_repeat = repeat;
        }
    }
}

} // namespace octopus

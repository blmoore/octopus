//
//  aligned_read.cpp
//  Octopus
//
//  Created by Daniel Cooke on 18/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "aligned_read.h"

AlignedRead::AlignedRead(const AlignedRead& other)
:
reference_region_ {other.reference_region_},
read_group_ {other.read_group_},
sequence_ {other.sequence_},
qualities_ {other.qualities_},
cigar_string_ {other.cigar_string_},
next_segment_ {((other.next_segment_ != nullptr) ?
                std::make_unique<NextSegment>(*other.next_segment_) : nullptr) },
flags_ {other.flags_},
mapping_quality_ {other.mapping_quality_},
is_compressed_ {other.is_compressed()}
{}

AlignedRead& AlignedRead::operator=(const AlignedRead& other)
{
    AlignedRead temp {other};
    swap(*this, temp);
    return *this;
}

void swap(AlignedRead& lhs, AlignedRead& rhs) noexcept
{
    std::swap(lhs.reference_region_, rhs.reference_region_);
    std::swap(lhs.read_group_, rhs.read_group_);
    std::swap(lhs.sequence_, rhs.sequence_);
    std::swap(lhs.cigar_string_, rhs.cigar_string_);
    std::swap(lhs.qualities_, rhs.qualities_);
    std::swap(lhs.next_segment_, rhs.next_segment_);
    std::swap(lhs.flags_, rhs.flags_);
    std::swap(lhs.mapping_quality_, rhs.mapping_quality_);
    std::swap(lhs.is_compressed_, rhs.is_compressed_);
}

//
// NextSegment public methods
//

const std::string& AlignedRead::NextSegment::get_contig_name() const
{
    return contig_name_;
}

AlignedRead::SizeType AlignedRead::NextSegment::get_begin() const noexcept
{
    return begin_;
}

AlignedRead::SizeType AlignedRead::NextSegment::get_inferred_template_length() const noexcept
{
    return inferred_template_length_;
}

bool AlignedRead::NextSegment::is_marked_unmapped() const
{
    return flags_[0];
}

bool AlignedRead::NextSegment::is_marked_reverse_mapped() const
{
    return flags_[1];
}

//
// AlignedRead public methods
//

const GenomicRegion& AlignedRead::get_region() const noexcept
{
    return reference_region_;
}

const AlignedRead::SequenceType& AlignedRead::get_sequence() const noexcept
{
    return sequence_;
}

const AlignedRead::Qualities& AlignedRead::get_qualities() const noexcept
{
    return qualities_;
}

void AlignedRead::zero_front_qualities(SizeType num_bases) noexcept
{
    std::for_each(std::begin(qualities_), std::begin(qualities_) + num_bases,
                  [] (auto& quality) { quality = 0; });
}

void AlignedRead::zero_back_qualities(SizeType num_bases) noexcept
{
    std::for_each(std::rbegin(qualities_), std::rbegin(qualities_) + num_bases,
                  [] (auto& quality) { quality = 0; });
}

AlignedRead::QualityType AlignedRead::get_mapping_quality() const noexcept
{
    return mapping_quality_;
}

AlignedRead::SizeType AlignedRead::get_sequence_size() const noexcept
{
    return static_cast<SizeType>(sequence_.size());
}

const CigarString& AlignedRead::get_cigar_string() const noexcept
{
    return cigar_string_;
}

const std::unique_ptr<AlignedRead::NextSegment>& AlignedRead::get_next_segment() const
{
    if (is_chimeric()) {
        return next_segment_;
    } else {
        throw std::runtime_error {"Read does not have a next segment"};
    }
}

AlignedRead::FlagData AlignedRead::get_flags() const
{
    FlagData flags {};
    flags.is_marked_multiple_read_template       = is_marked_multiple_read_template();
    flags.is_marked_all_segments_in_read_aligned = is_marked_all_segments_in_read_aligned();
    flags.is_marked_unmapped                     = is_marked_unmapped();
    flags.is_marked_reverse_mapped               = is_marked_reverse_mapped();
    flags.is_marked_secondary_alignment          = is_marked_secondary_alignment();
    flags.is_marked_qc_fail                      = is_marked_qc_fail();
    flags.is_marked_duplicate                    = is_marked_duplicate();
    flags.is_marked_supplementary_alignment      = is_marked_supplementary_alignment();
    
    return flags;
}

bool AlignedRead::is_chimeric() const noexcept
{
    return next_segment_ != nullptr;
}

bool AlignedRead::is_marked_all_segments_in_read_aligned() const
{
    return flags_[0];
}

bool AlignedRead::is_marked_multiple_read_template() const
{
    return flags_[1];
}

bool AlignedRead::is_marked_unmapped() const
{
    return flags_[2];
}

bool AlignedRead::is_marked_reverse_mapped() const
{
    return flags_[3];
}

bool AlignedRead::is_marked_secondary_alignment() const
{
    return flags_[4];
}

bool AlignedRead::is_marked_qc_fail() const
{
    return flags_[5];
}

bool AlignedRead::is_marked_duplicate() const
{
    return flags_[6];
}

bool AlignedRead::is_marked_supplementary_alignment() const
{
    return flags_[7];
}

//
// Private methods
//

bool AlignedRead::is_compressed() const noexcept
{
    return is_compressed_;
}

void AlignedRead::set_compressed() noexcept
{
    is_compressed_ = true;
}

void AlignedRead::set_uncompressed() noexcept
{
    is_compressed_ = false;
}

AlignedRead::Flags AlignedRead::get_flags(const FlagData& flags)
{
    Flags result {};
    result[0] = flags.is_marked_all_segments_in_read_aligned;
    result[1] = flags.is_marked_multiple_read_template;
    result[2] = flags.is_marked_unmapped;
    result[3] = flags.is_marked_reverse_mapped;
    result[4] = flags.is_marked_secondary_alignment;
    result[5] = flags.is_marked_qc_fail;
    result[6] = flags.is_marked_duplicate;
    result[7] = flags.is_marked_supplementary_alignment;
    return result;
}

AlignedRead::NextSegment::Flags AlignedRead::NextSegment::get_flags(const FlagData& flags)
{
    Flags result {};
    result[0] = flags.is_marked_unmapped;
    result[1] = flags.is_marked_reverse_mapped;
    return result;
}

// Non-member methods

AlignedRead splice(const AlignedRead& read, const GenomicRegion& region)
{
    if (!contains(read, region)) {
        throw std::runtime_error {"cannot splice AlignedRead region that is not contained"};
    }
    
    auto reference_offset = get_begin(region) - get_begin(read);
    
    auto uncontained_cigar_splice = reference_splice(read.get_cigar_string(), GenomicRegion::SizeType {}, reference_offset);
    auto contained_cigar_splice   = reference_splice(read.get_cigar_string(), reference_offset, size(region));
    
    auto sequence_offset = sequence_size(uncontained_cigar_splice);
    auto sequence_length = sequence_size(contained_cigar_splice);
    
    AlignedRead::SequenceType sequence_splice(read.get_sequence().cbegin() + sequence_offset,
                                              read.get_sequence().cbegin() + sequence_offset + sequence_length);
    AlignedRead::Qualities qualities_splice(read.get_qualities().cbegin() + sequence_offset,
                                            read.get_qualities().cbegin() + sequence_offset + sequence_length);
    
    return AlignedRead {
        region,
        std::move(sequence_splice),
        std::move(qualities_splice),
        std::move(contained_cigar_splice),
        read.get_mapping_quality(),
        read.get_flags()
    };
}
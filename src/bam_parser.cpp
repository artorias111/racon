/*!
 * @file bam_parser.cpp
 *
 * @brief BAM parser source file using htslib
 */

#include <cstring>
#include <algorithm>

#include "bam_parser.hpp"
#include "overlap.hpp"

namespace racon {

BamParser<Overlap>::BamParser(const std::string& path)
        : path_(path), file_(nullptr), header_(nullptr), alignment_(nullptr),
          num_objects_read_(0), is_eof_(false) {

    file_ = hts_open(path_.c_str(), "r");
    if (file_ == nullptr) {
        fprintf(stderr, "[racon::BamParser] error: unable to open file %s!\n",
            path_.c_str());
        exit(1);
    }

    header_ = sam_hdr_read(file_);
    if (header_ == nullptr) {
        fprintf(stderr, "[racon::BamParser] error: unable to read header from %s!\n",
            path_.c_str());
        hts_close(file_);
        exit(1);
    }

    alignment_ = bam_init1();
    if (alignment_ == nullptr) {
        fprintf(stderr, "[racon::BamParser] error: unable to allocate alignment!\n");
        bam_hdr_destroy(header_);
        hts_close(file_);
        exit(1);
    }
}

BamParser<Overlap>::~BamParser() {
    if (alignment_ != nullptr) {
        bam_destroy1(alignment_);
    }
    if (header_ != nullptr) {
        bam_hdr_destroy(header_);
    }
    if (file_ != nullptr) {
        hts_close(file_);
    }
}

void BamParser<Overlap>::reset() {
    if (file_ != nullptr) {
        hts_close(file_);
    }

    file_ = hts_open(path_.c_str(), "r");
    if (file_ == nullptr) {
        fprintf(stderr, "[racon::BamParser] error: unable to reopen file %s!\n",
            path_.c_str());
        exit(1);
    }

    if (header_ != nullptr) {
        bam_hdr_destroy(header_);
    }

    header_ = sam_hdr_read(file_);
    if (header_ == nullptr) {
        fprintf(stderr, "[racon::BamParser] error: unable to read header from %s!\n",
            path_.c_str());
        hts_close(file_);
        exit(1);
    }

    num_objects_read_ = 0;
    is_eof_ = false;
}

bool BamParser<Overlap>::parse(std::vector<std::unique_ptr<Overlap>>& dst,
    uint64_t max_bytes) {

    if (is_eof_) {
        return false;
    }

    uint64_t current_bytes = 0;
    uint64_t num_objects = 0;

    while (sam_read1(file_, header_, alignment_) >= 0) {
        uint32_t flag = alignment_->core.flag;

        // Skip unmapped reads and secondary/supplementary alignments
        if ((flag & BAM_FUNMAP) || (flag & BAM_FSECONDARY) || (flag & BAM_FSUPPLEMENTARY)) {
            continue;
        }

        // Get query name
        const char* q_name = bam_get_qname(alignment_);
        uint32_t q_name_length = alignment_->core.l_qname - 1;

        // Get target name
        const char* t_name = header_->target_name[alignment_->core.tid];
        uint32_t t_name_length = strlen(t_name);

        // Get target position (1-based in BAM, convert to 1-based for SAM constructor)
        uint32_t t_begin = alignment_->core.pos + 1;

        // Get mapping quality
        uint32_t mapping_quality = alignment_->core.qual;

        // Build CIGAR string
        std::string cigar_str;
        uint32_t* cigar = bam_get_cigar(alignment_);
        for (uint32_t i = 0; i < alignment_->core.n_cigar; ++i) {
            cigar_str += std::to_string(bam_cigar_oplen(cigar[i]));
            cigar_str += bam_cigar_opchr(cigar[i]);
        }

        // Get sequence (only needed for length calculation, not stored)
        uint32_t sequence_length = alignment_->core.l_qseq;

        // Get quality (optional, may be empty)
        uint8_t* qual = bam_get_qual(alignment_);
        std::string quality_str;
        if (qual[0] != 0xff) {
            quality_str.resize(sequence_length);
            for (uint32_t i = 0; i < sequence_length; ++i) {
                quality_str[i] = static_cast<char>(qual[i] + 33);
            }
        }

        // Create Overlap using the SAM constructor
        // The constructor expects: q_name, q_name_length, flag, t_name, t_name_length,
        // t_begin, mapping_quality, cigar, cigar_length, t_next_name, t_next_name_length,
        // t_next_begin, template_length, sequence, sequence_length, quality, quality_length
        dst.emplace_back(new Overlap(
            q_name, q_name_length,
            flag,
            t_name, t_name_length,
            t_begin,
            mapping_quality,
            cigar_str.c_str(), static_cast<uint32_t>(cigar_str.length()),
            "*", 1,  // t_next_name (not used)
            0,       // t_next_begin (not used)
            0,       // template_length (not used)
            nullptr, sequence_length,  // sequence data not needed, just length
            quality_str.empty() ? nullptr : quality_str.c_str(),
            static_cast<uint32_t>(quality_str.length())
        ));

        ++num_objects;
        ++num_objects_read_;

        // Estimate bytes read (approximate)
        current_bytes += q_name_length + t_name_length + cigar_str.length() +
                        sequence_length + quality_str.length() + 100;

        if (current_bytes >= max_bytes) {
            return true;
        }
    }

    is_eof_ = true;
    return num_objects > 0;
}

std::unique_ptr<BamParser<Overlap>> createBamParser(const std::string& path) {
    return std::unique_ptr<BamParser<Overlap>>(new BamParser<Overlap>(path));
}

}

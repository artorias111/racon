/*!
 * @file bam_parser.hpp
 *
 * @brief BAM parser header file using htslib
 */

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <atomic>

#include <htslib/sam.h>
#include <htslib/hts.h>

namespace racon {

template<class T>
class BamParser;

class Overlap;

template<>
class BamParser<Overlap> {
public:
    BamParser(const std::string& path);
    BamParser(const BamParser&) = delete;
    BamParser& operator=(const BamParser&) = delete;
    ~BamParser();

    void reset();

    bool parse(std::vector<std::unique_ptr<Overlap>>& dst, uint64_t max_bytes);

private:
    std::string path_;
    htsFile* file_;
    bam_hdr_t* header_;
    bam1_t* alignment_;
    std::atomic<uint64_t> num_objects_read_;
    bool is_eof_;
};

std::unique_ptr<BamParser<Overlap>> createBamParser(const std::string& path);

}

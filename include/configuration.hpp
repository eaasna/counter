// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <filesystem>
#include <seqan3/search/kmer_index/shape.hpp>

struct configuration
{
    std::filesystem::path fasta_input{};
    std::filesystem::path output{"output_index"};
    bool verbose{}; // Default is false.
    std::string shape_input{"1111"};
    uint8_t shape_size;
    seqan3::shape shape;
    uint8_t window_input{8};

};

// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <gtest/gtest.h>

#include <get_kmer_from_hash.hpp>

#include "app_test.hpp"

// To prevent issues when running multiple API tests in parallel, give each API test unique names:
struct hash_to_kmer : public app_test
{};

TEST_F(hash_to_kmer, ungapped)
{
    uint64_t hash{28};
    uint8_t kmer_len{3};

    std::string expected{"CTA"};
    seqan3::shape shape{seqan3::ungapped{kmer_len}};
    auto actual = get_kmer_from_hash(shape, hash);

    EXPECT_EQ(expected, actual);
}

TEST_F(hash_to_kmer, lowest)
{
    uint64_t hash{0};
    uint8_t kmer_len{3};

    std::string expected{"AAA"};
    seqan3::shape shape{seqan3::ungapped{kmer_len}};
    auto actual = get_kmer_from_hash(shape, hash);

    EXPECT_EQ(expected, actual);
}

TEST_F(hash_to_kmer, highest)
{
    uint64_t hash{63};
    uint8_t kmer_len{3};

    std::string expected{"TTT"};
    seqan3::shape shape{seqan3::ungapped{kmer_len}};
    auto actual = get_kmer_from_hash(shape, hash);

    EXPECT_EQ(expected, actual);
}

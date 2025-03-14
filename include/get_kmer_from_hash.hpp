// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <seqan3/search/kmer_index/shape.hpp>

std::string get_kmer_from_hash(seqan3::shape const & shape, uint64_t const hash);

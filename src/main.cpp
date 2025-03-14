// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "configuration.hpp"
#include "get_shape.hpp"
#include "counting_index.hpp"

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <sharg/all.hpp>

#include <cereal/archives/xml.hpp>
#include <iostream>

using namespace seqan3::literals;

int main(int argc, char ** argv)
{
    // Configuration
    configuration config{};

    // Parser
    sharg::parser parser{"K-mer-Counter", argc, argv};

    // General information.
    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    // Positional option: The FASTQ file to convert.
    parser.add_positional_option(config.fasta_input,
                                 sharg::config{.description = "The FASTA file with sequences to count.",
                                               .validator = sharg::input_file_validator{{"fa", "fasta"}}});
    
    
    // Output file. Default: print to terminal 
    parser.add_option(config.output,
                      sharg::config{.short_id = 'o',
                                    .long_id = "output",
                                    .description = "The output file with counted Kmere.",
                                    .default_message = "Print to terminal (stdout)",
                                    .validator = sharg::output_file_validator{}});

    //Validator for Kmer shape
    sharg::regex_validator shape_validator{"[0-1]"};

    //Kmer shape
    parser.add_option(config.shape_input,
                        sharg::config{.short_id = 's', 
                                      .long_id = "shape", 
                                      .description = "The shape of the Kmer.",
                                      .validator = shape_validator});
    


    //Window size input
    parser.add_option(config.window_input, 
                        sharg::config{.short_id = 'w',
                                      .long_id = "window", 
                                      .description = "The window size for the Kmer"});

                                   
    // Flag: Verose output.
    parser.add_flag(
        config.verbose,
        sharg::config{.short_id = 'v', .long_id = "verbose", .description = "Give more detailed information."});
    

    try
    {
        parser.parse(); // Trigger command line parsing.
    }
    catch (sharg::parser_error const & ext) // Catch user errors.
    {
        std::cerr << "Parsing error. " << ext.what() << '\n'; // Give error message.
        return -1;
    }


    //auto sequence_size  ...
    config.shape = get_shape(config.shape_input);
    config.shape_size = config.shape.size();

    //Validating Window Input
    if(config.window_input <= config.shape_size){
        throw std::runtime_error("Window size to big");

    }

    counting_index kmer_index(config);
    std::ofstream os(config.output);
    cereal::XMLOutputArchive archive( os );
    kmer_index.save(archive);

    if (config.verbose) // If flag is set.
    	std::cerr << "Counting was a success. Congrats!\n";

    return 0;
}

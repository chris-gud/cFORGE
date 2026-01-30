// Primer.hpp
#ifndef PRIMER_H
#define PRIMER_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <cmath>
#include <cctype>
#include <stdexcept>
#include "nucleotide_utils.hpp"

struct Primer {

    std::string seq;                        // Primer sequence. MUST be 5' -> 3'
    size_t length;                          // Length of sequence
    unsigned long wallace_melting_temp;     // Melting temp using Wallace approximation. Units: °C
    double santalucia_melting_temp;         // Melting temp using SantaLucia 2004. Units: °C
    size_t num_A = 0;                       // Number of Adenines in sequence
    size_t num_T = 0;                       // Number of Thymines in sequence
    size_t num_G = 0;                       // Number of Guanines in sequence
    size_t num_C = 0;                       // Number of Cytosines in sequence
    float GC_content;                       // Proportion of G's and C's in sequence. GC_content ∈ [0.0 , 1.0]
    bool GC_clamp;                          // True if exists a G or C within 5 nucleotides of 3' end. False if not
    bool homo_dimerizes;                    // True if dimerizes_with another primer of the same sequence
    size_t max_length_homopolymer_run;      // The length of the longest homopolymer run ("AAAAA") in the sequence
    bool hairpin;                           // True if primer is likely to form deleterious hairpin structure(s)
    
    /*
     * Primer constructor (std::string sequence)
     */
    Primer (std::string sequence);
    
    
    /*
     * Returns true if Primer this and Primer other are LIKELY to dimerize.
     * Good for reuse to determine homodimerization OR heterodimerization.
     * Utilizes the dimerize() function in nucleotide_utils.
     * returns true if dimerize(this->seq, other_seq) returns a double >= 0.75
     */
    bool dimerizes_with(const std::string& other_seq);
    
    
    /*
     * bool is_good() returns true if this Primer satisfies the (weak) primer design rules as notated below:
     *
     * Primer rules:
     *   1.) 18-25 bases long ideally. PCR mutagenesis allows up to 45 n.t.
     *   2.) 40-70% GC content
     *   3.) Melting temperature between 55-65°C
     *   4.) GC clamp at the 3' end
     *   5.) Avoid primer homodimerization, homopolymer runs, hairpins, off-target binding, heterodimerization
     */
    bool is_good();
    
    
    /*
     * bool binds_off_target(std::string &gene_seq) returns true if the reverse complement of the
     *  primer sequence is found exactly in gene_seq OR the reverse complement of the gene_seq MORE THAN ONCE.
     * For use within main.
     * SHOULD BE IMPROVED to find imperfect matches.
     */
    bool binds_off_target(const std::string &gene_seq);
    
    
    /*
     * to_string() returns information about Primer in std::string format
     */
    std::string to_string();
    
private:
    /*
     * bool forms_hairpin returns true if Primer sequence is likely to form a deleterious secondary structure.
     * Utilizes the hairpin() function in nucleotide_utils.
     * returns true if hairpin(this->seq) returns a double >= 0.65
     */
    bool forms_hairpin();
    
    /* 
     * Calculates the melting temperature of the Primer using the SantaLucia (2004) algorithm.
     */
    double calculate_tm_santalucia_2004(
        const std::string& primer,
        double primer_conc = 5e-7,  // 500 nM
        double na_conc = 0.05       // 50 mM Na+
    );

};
#endif

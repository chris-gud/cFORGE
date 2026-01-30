// nucleotide_utils.hpp
#pragma once

#include <string>
#include <cctype>
#include <unordered_map>
#include <algorithm>
#include <cmath>

namespace nucleotide {

    /*
     * Normalizes a nucleic acid sequence; removes spaces, non-standard characters, and capitalizes characters.
     * Supports nucleotides A,T,G,C,U,N
     */
    std::string normalize(const std::string& seq);

    /*
     * Returns the char complement of a nucleotide base (A↔T, C↔G).
     * Supports lowercase nucleotides (insertion/substitution/duplication mutated) and char '-' (deleted nucleotide).
     * Returns '?' for undefined characters.
     */
    char complement_of(char base);

    /*
     * Returns the complement of a std::string nucleotide sequence.
     */
    std::string complement_of(const std::string &sequence);

    /*
     * Returns the reverse complement of a std::string nucleotide sequence.
     */
    std::string revcomp_of(const std::string &sequence);

    /*
     * Replaces the "==" operator use case when seeing if a character in a main_sequence satisfies a non-standard
     *  nucleic acid code.
     * For support of non-standard nucleic acid codes in char rec_seq_char
     *
     * Ex. 1: matches(rec_seq_char = 'K', main_seq_char = 'G')
     *       returns TRUE since G is a pyrimidine (K)
     * Ex. 2: matches(rec_seq_char = 'M', main_seq_char = 'T')
     *       returns FALSE since T is not an amino group nucleotide (M)
     */
    bool matches(char rec_seq_char, char main_seq_char);

    /*
     * New find method for use in nucleic acid recognition sequence detection.
     * Supports non-standard nucleic acid codes:
     *  N = any base (A, T, G, C)
     *  Y = pyrimidine (C or T)
     *  R = purine (A or G)
     *  M = amino (A or C) (NH2 group)
     *  K = keto (G or T) (C=O group)
     *
     * Ex. 1: find_DNA_seq(recognition_seq = "AR", main_seq = "ACTAA")
     *        returns 3 since the substring "AA" is recognized as satisfying the recognition sequence "AR"
     */
    size_t find_DNA_seq(const std::string& recognition_seq, const std::string& main_seq);

    /*
     * Function for mutating a gene sequence given a notated mutation.
     * Returns the mutated gene sequence.
     * Mutated nucleotides/sequences are *lowercase* in returned mutated sequence.
     * Deleted nucleotides/sequences are indicated by one "-" in returned mutated sequence.
     * param gene_seq remains unchanged.
     * param mutation is given in mutation notation as notated by https://atlasgeneticsoncology.org/teaching/30067/nomenclature-for-the-description-of-mutations-and-other-sequence-variations
     */
    std::string mutate(const std::string& gene_seq, const std::string& mutation);

    /*
     * extract_old_sequence() returns the old sequence to be mutated given a mutation in standard notation.
     * param gene_seq remains unchanged.
     * param mutation: see function comment for nucleotide::mutate())
     * If:
     *      substitution mutation: returns the old nucleotide to be sustituted. (ex. "c.123A>G" --> returns "A")
     *      deletion mutation: returns old sequence that was there "". (ex. "c.546delT" --> returns "T")
     *      duplication mutation: returns the sequence to be duplicated. (ex. "c.546dupT" --> returns "T")
     *      insertion mutation: returns empty string. (ex. "c.186_187insGCGTGA" --> returns "")
     */
    std::string extract_old_sequence(const std::string& gene_seq, const std::string& mutation);

    /*
     * extract_mutation_sequence() returns the new mutated sequence given a mutation in standard notation.
     * It does NOT return the FULL mutated GENE sequence. For that function, see nucleotide::mutate().
     * param gene_seq remains unchanged.
     * param mutation: see function comment for nucleotide::mutate())
     * If:
     *      substitution mutation: returns the new substituted nucleotide. (ex. "c.123A>G" --> returns "G")
     *      deletion mutation: returns empty string "". (ex. "c.546delT" --> returns "")
     *      duplication mutation: returns the sequence, duplicated. (ex. "c.546dupT" --> returns "TT")
     *      insertion mutation: returns the sequence to be inserted. (ex. "c.186_187insGCGTGA" --> returns "GCGTGA")
     */
    std::string extract_mutation_sequence(const std::string& gene_seq, const std::string& mutation);

    /*
     * Determines and returns the 0-BASED index of the FIRST nucleotide edited by a given mutation.
     * param mutation is given in standard mutation notation (see function comment nucleotide::mutate()).
     * Ex. extract_mutation_position("c.586_591dup") returns 567 (0-based indexing subtracts 1 from 1-based)
     */
    size_t extract_mutation_position(const std::string& mutation);

    /*
     * Estimates the "likelihood" (0.0 to 1.0) that two nucleic acid sequences will form a stable
     * intermolecular dimer (e.g., primer–primer dimerization).
     *
     * The function aligns one sequence against the reverse complement of the other
     * and scans all relative offsets to identify contiguous Watson–Crick complementary
     * regions. Candidate duplexes are evaluated using a simplified nearest-neighbor
     * thermodynamic model, in which base-stacking interactions determine the relative
     * Gibbs free energy (ΔG) of duplex formation.
     *
     * Nearest-neighbor parameters and the thermodynamic framework are based on:
     *   SantaLucia, J. (1998) "A unified view of polymer, dumbbell, and oligonucleotide
     *   DNA nearest-neighbor thermodynamics." Proc. Natl. Acad. Sci. USA.
     *   https://www.pnas.org/doi/10.1073/pnas.95.4.1460
     *
     * This algorithm places additional emphasis on complementarity near the 3′ ends
     * of both sequences, as 3′-end pairing can be extended by DNA polymerase and is a
     * primary cause of primer-dimer artifacts in PCR. This design choice mirrors the
     * screening logic used in common primer design and analysis tools, including:
     *   Primer3: https://primer3.org/
     *   Untergasser et al. (2012) Nucleic Acids Res.
     *   https://academic.oup.com/nar/article/40/15/e115/2412948
     *
     * Practical ΔG-based thresholds and GC-rich penalties are consistent with guidance
     * from oligonucleotide analysis tools such as IDT OligoAnalyzer:
     *   https://www.idtdna.com/pages/tools/oligoanalyzer
     *
     * The resulting features (minimum ΔG across candidate duplexes, longest contiguous
     * complementary segment, and 3′-end complementarity metrics) are combined into a
     * heuristic likelihood value in the range [0.0, 1.0]. Larger values indicate a
     * higher relative risk of biologically relevant dimerization. The returned value
     * is a comparative risk estimate rather than an absolute thermodynamic probability.
     */
    double dimerize(std::string seq1, std::string seq2);

    /*
     * Estimates the "likelihood" [0.0, 1.0] that a single nucleic acid sequence will form a stable
     * intramolecular hairpin (stem–loop) structure.
     *
     * Method summary:
     *   - Normalizes to DNA alphabet (uppercases; U->T). Non-standard bases are unsupported and will
     *     break potential stems.
     *   - Enumerates plausible hairpin geometries by selecting a loop length and a left-stem start,
     *     then searches for the longest contiguous reverse-complementary "stem" across the loop.
     *   - Scores candidate stems using a simplified SantaLucia-style nearest-neighbor ΔG°37 stack sum
     *     (no salt/concentration/initiation corrections). Longer and GC-rich stems are treated as
     *     more stable. Stems involving the 3′ end are weighted more heavily (PCR relevance).
     *
     * Returns:
     *   - A heuristic likelihood in [0.0, 1.0] (higher = higher relative hairpin risk).
     *     NOTE: This is a comparative screening score, not an absolute physical probability.
     */
    double hairpin(std::string seq);


} // namespace nucleotide

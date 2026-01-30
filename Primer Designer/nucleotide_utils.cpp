// nucleotide_utils.cpp
#include "nucleotide_utils.hpp"

namespace nucleotide {
    // Anonymous namespace for functions for use ONLY within namespace nucleotide
    namespace {
        /*
         * extract_num_from_str extracts the first contiguous number from a string.
         * For INTERNAL use for determining the number within a mutation string.
         * Example 1: extract_location(s = "c.123A>G") returns 123.
         * Example 2: extract_location(s = "c.1086_1087insGCGTGA") returns 1086.
         * Param s MUST contain ONLY ONE contiguous number.
         */
        size_t extract_num_from_str(const std::string &s) {
            size_t value = 0;
            bool reading = false;

            for (char c : s) {
                if (std::isdigit(static_cast<unsigned char>(c))) {
                    reading = true;
                    value = value * 10 + (c - '0');
                } else if (reading) {
                    break;
                }
            }
            return value;
        }
    
        /*
         * to_lower_str makes a string lowercase.
         * For internal use only (used in function mutate).
         */
        std::string to_lower_str(const std::string &s) {
            std::string lower = "";
            for (char c : s) {
                lower += (char)std::tolower(c);
            }
            return lower;
        }
    } // End anonymous namespace


    /*
     * Normalizes a nucleic acid sequence; removes spaces, non-standard characters, and capitalizes characters.
     * Supports nucleotides A,T,G,C,U,N
     */
    std::string normalize(const std::string& seq) {
        std::string out;
        out.reserve(seq.size());

        for (char c : seq) {
            if (std::isspace(static_cast<unsigned char>(c)))
                continue;

            char base = std::toupper(static_cast<unsigned char>(c));
            switch (base) {
                case 'A':
                case 'C':
                case 'G':
                case 'T':
                case 'U':
                case 'N':
                    out.push_back(base);
                    break;
                default:
                    // drop non-standard characters
                    break;
            }
        }
        return out;
    }


    /*
     * Returns the char complement of a nucleotide base (A↔T, C↔G).
     * Supports lowercase nucleotides (insertion/substitution/duplication mutated) and char '-' (deleted nucleotide).
     * Returns '?' for undefined characters.
     */
    char complement_of(char c) {
        switch (c) {
            case 'A': return 'T';
            case 'T': return 'A';
            case 'G': return 'C';
            case 'C': return 'G';
            case '-': return '-';
            case 'a': return 't';
            case 't': return 'a';
            case 'g': return 'c';
            case 'c': return 'g';
            default:  return '?';
        }
    }

    /*
     * Returns the complement of a std::string nucleotide sequence.
     */
    std::string complement_of(const std::string& sequence) {
        std::string comp_seq;
        comp_seq.reserve(sequence.size());

        for (char c : sequence) {
            comp_seq.push_back(complement_of(c));
        }

        return comp_seq;
    }

    /*
     * Returns the reverse complement of a std::string nucleotide sequence.
     */
    std::string revcomp_of(const std::string& sequence) {
        std::string revcomp;
        revcomp.reserve(sequence.size());

        for (size_t i = sequence.size(); i-- > 0; ) {
            revcomp.push_back(complement_of(sequence[i]));
        }

        return revcomp;
    }

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
    bool matches(char rec_seq_char, char main_seq_char) {
        switch (rec_seq_char) {
            case 'A': return main_seq_char == 'A';
            case 'T': return main_seq_char == 'T';
            case 'G': return main_seq_char == 'G';
            case 'C': return main_seq_char == 'C';
            case 'N': return true;
            case 'Y': return main_seq_char == 'C' || main_seq_char == 'T';
            case 'R': return main_seq_char == 'A' || main_seq_char == 'G';
            case 'M': return main_seq_char == 'A' || main_seq_char == 'C';
            case 'K': return main_seq_char == 'G' || main_seq_char == 'T';
            case 'S': return main_seq_char == 'G' || main_seq_char == 'C'; // strong
            case 'W': return main_seq_char == 'A' || main_seq_char == 'T'; // weak
            case 'B': return main_seq_char == 'C' || main_seq_char == 'G' || main_seq_char == 'T'; // not A
            case 'D': return main_seq_char == 'A' || main_seq_char == 'G' || main_seq_char == 'T'; // not C
            case 'H': return main_seq_char == 'A' || main_seq_char == 'C' || main_seq_char == 'T'; // not G
            case 'V': return main_seq_char == 'A' || main_seq_char == 'C' || main_seq_char == 'G'; // not T
            default:  return false;
        }
    }

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
    size_t find_DNA_seq(const std::string& recognition_seq, const std::string& main_seq){
        if (recognition_seq.size() > main_seq.size())
            return std::string::npos;

        for (size_t i = 0; i <= main_seq.size() - recognition_seq.size(); i++) {
            bool match = true;

            for (size_t k = 0; k < recognition_seq.size(); k++) {
                if (!matches(recognition_seq[k], main_seq[i + k])) {
                    match = false;
                    break;
                }
            }
            
            if (match)
                return i;
        }

        return std::string::npos;
    }


    /*
     * Function for mutating a gene sequence given a notated mutation.
     * Returns the mutated gene sequence.
     * Mutated nucleotides/sequences are *lowercase* in returned mutated sequence.
     * Deleted nucleotides/sequences are indicated by one "-" in returned mutated sequence.
     * param gene_seq remains unchanged.
     * param mutation is given in mutation notation as notated by https://atlasgeneticsoncology.org/teaching/30067/nomenclature-for-the-description-of-mutations-and-other-sequence-variations
     */
    std::string mutate(const std::string& gene_seq, const std::string& mutation){
        std::string normalized_gene_seq = normalize(gene_seq);
        std::string mutated_gene_seq = "";
        
        // Determine reference sequence.
        bool coding = mutation.find('c') != std::string::npos;
        bool genomic = mutation.find('g') != std::string::npos;
        bool mitochondrial = mutation.find('m') != std::string::npos;
        bool protein = mutation.find('p') != std::string::npos;
        bool rna = mutation.find('r') != std::string::npos;
        
        if (coding || genomic || mitochondrial) {   // All are DNA reference sequence mutations.
            
            // Determine mutation type.
            bool substitution = mutation.find('>') != std::string::npos;
            bool deletion = mutation.find("del") != std::string::npos;
            bool duplication = mutation.find("dup") != std::string::npos;
            bool insertion = mutation.find("ins") != std::string::npos;
            
            // Perform given mutation.
            if (substitution) {
                /* Format of substitution mutation:
                    "c.123A>G"
                 */
                
                // Find substitution index. In example above: 122
                // Subtract 1 to convert to 0-based indexing.
                std::size_t position_mutation = extract_num_from_str(mutation) - 1; // 0-BASED
                
                // Make sure given index is not > than gene_seq size.
                if (position_mutation > normalized_gene_seq.size()) {
                    std::string error_message = "";
                    error_message += "?Error: Insertion mutation index is greater than the length of the gene sequence.";
                    return error_message;
                }
                
                /* Determine original nucleotide and mutated mucleotide.
                    original_nucleotide will be to the left of '>', mutated_nucleotide will be to the right.
                 */
                // Index of '>' within string mutation
                size_t index_greaterthan = mutation.find('>');
                size_t index_original_nucleotide = index_greaterthan - 1;
                size_t index_mutated_nucleotide = index_greaterthan + 1;
                char original_nucleotide = mutation[index_original_nucleotide];
                char mutated_nucleotide = mutation[index_mutated_nucleotide];
                
                /* Generate substitution mutated gene sequence. */
                mutated_gene_seq = normalized_gene_seq;
                char nucleotide_at_pos = mutated_gene_seq[position_mutation];
                if (nucleotide_at_pos == original_nucleotide) {
                    // Mutate.
                    mutated_gene_seq[position_mutation] = (char)std::tolower(mutated_nucleotide);
                    // Substitution mutation success.
                    return mutated_gene_seq;
                } else {
                    // return mutation notation incorrect error
                    size_t position_mutation_1based = position_mutation + 1;
                    std::string error_message = "";
                    error_message += "?Error: mutation notation incorrect. Nucleotide at position ";
                    error_message += std::to_string(position_mutation_1based) + " is " + nucleotide_at_pos;
                    error_message += std::string(", not ") + original_nucleotide + std::string(".");
                    return error_message;
                }
                
            } else if (deletion) {
                /* Formats of deletion mutations: 
                    "c.546delT"  --> single base at position 546 deleted.
                    "c.586_588del"  --> bases at positions 586, 587, 588 deleted.
                 */
                
                // Determine if single base deletion or multiple.
                bool multiple_deletions = mutation.find('_') != std::string::npos;
                if (multiple_deletions) {   // Multiple base deletion
                    /* Deal with indexing. */
                    // Generate substring of mutation to include only FIRST number.
                    std::string mutation_substr_1 = mutation.substr(0, mutation.find('_'));
                    // Generate substring of mutation to include only SECOND number.
                    std::string mutation_substr_2 = mutation.substr(mutation.find('_'));
                    // Find index of first nucleotide to be deleted.
                    size_t position_first = extract_num_from_str(mutation_substr_1) - 1; //0-BASED
                    // Find index of last nucleotide to be deleted.
                    size_t position_last = extract_num_from_str(mutation_substr_2) - 1; //0-BASED
                    // Calculate number of nucleotides to be deleted.
                    size_t num_deletions = position_last - position_first + 1;
                    
                    // Make sure given index is not > than gene_seq size.
                    if (position_last > normalized_gene_seq.size()) {
                        std::string error_message = "";
                        error_message += "?Error: Deletion mutation index is greater than the length of the gene sequence.";
                        return error_message;
                    }
                    
                    /* Generate multiple deletion mutated gene sequence. */
                    mutated_gene_seq = normalized_gene_seq;
                    // Delete num_deletions chars at index position_first.
                    mutated_gene_seq.erase(position_first, num_deletions);
                    // Add "-" to location of deleted nucleotides.
                    mutated_gene_seq = mutated_gene_seq.substr(0,position_first) + "-" + mutated_gene_seq.substr(position_first);
                    
                    // Multiple deletions mutation success.
                    return mutated_gene_seq;
                    
                } else {    // Single base deletion.
                    // Find deletion index.
                    // Subtract 1 to convert to 0-based indexing.
                    size_t position_mutation = extract_num_from_str(mutation) - 1; // 0-BASED
                    
                    // Make sure given index is not > than gene_seq size.
                    if (position_mutation > normalized_gene_seq.size()) {
                        std::string error_message = "?Error: Deletion mutation index is greater than the length of the gene sequence.";
                        return error_message;
                    }
                    
                    // Generate single deletion mutated gene sequence.
                    mutated_gene_seq = normalized_gene_seq;
                    mutated_gene_seq.erase(position_mutation, 1); // Delete one char at index mutated_gene_seq.
                    // Add "-" to location of deleted nucleotide.
                    mutated_gene_seq = mutated_gene_seq.substr(0,position_mutation) + "-" + mutated_gene_seq.substr(position_mutation);
                    
                    // Single deletion mutation success.
                    return mutated_gene_seq;
                }
                
            } else if (duplication) {
                /* Formats of duplication mutation:
                    "c.546dupT"  --> duplicate the nucleotide at position 546.
                    "c.586_588dup"  --> duplicate the nucleotide sequence 586, 587, 588.
                 */
                
                // Determine if single base or sequence duplication.
                bool sequence_duplication = mutation.find('_') != std::string::npos;
                
                /* Identify variables. */
                size_t position_first;
                size_t position_last;
                size_t length_duplication_sequence;
                
                if (!sequence_duplication) { // Single base duplication.
                    /* Deal with indexing for single base duplication.*/
                    // Find index of nucleotide in sequence-to-be-duplicated.
                    position_first = extract_num_from_str(mutation) - 1; // 0-BASED
                    position_last = position_first;
                    
                    // 1 nucleotide to be duplicated.
                    length_duplication_sequence = 1;
                    
                } else { // Sequence duplication.
                    /* Deal with indexing for sequence duplication.*/
                    // Generate substring of mutation to include only FIRST number.
                    std::string mutation_substr_1 = mutation.substr(0, mutation.find('_'));
                    // Generate substring of mutation to include only SECOND number.
                    std::string mutation_substr_2 = mutation.substr(mutation.find('_'));
                    // Find index of first nucleotide in sequence-to-be-duplicated.
                    position_first = extract_num_from_str(mutation_substr_1) - 1; // 0-BASED
                    // Find index of last nucleotide in sequence-to-be-duplicated.
                    position_last = extract_num_from_str(mutation_substr_2) - 1; // 0-BASED
                    // Calculate number of nucleotides to be duplicated.
                    length_duplication_sequence = position_last - position_first + 1;
                }
                
                // Make sure given index is not > than gene_seq size.
                if (position_last > normalized_gene_seq.size()) {
                    std::string error_message = "";
                    error_message += "?Error: Duplication mutation index is greater than the length of the gene sequence.";
                    return error_message;
                }
                
                // Generate sequence-to-be-duplicated.
                std::string sequence_to_be_duplicated = normalized_gene_seq.substr(position_first, length_duplication_sequence);
                
                /* Generate duplication mutated gene sequence. */
                mutated_gene_seq = normalized_gene_seq.substr(0, position_first);
                mutated_gene_seq += sequence_to_be_duplicated;
                mutated_gene_seq += to_lower_str(sequence_to_be_duplicated); /* newly duplicated sequence should be lowercase */
                mutated_gene_seq += normalized_gene_seq.substr(position_last + 1);

                // Duplication mutation successful.
                return mutated_gene_seq;
                
            } else if (insertion) {
                /* Format of insertion mutation:
                    "c.546_547insT"  --> insertion of T between 546 and 547
                    "c.186_187insGCGTGA"  --> insertion of GCGTGA between 186 and 187
                 */
                
                // Find index of insertion.
                size_t position_insertion = extract_mutation_position(mutation); // 0-BASED.
                
                // Make sure given index is not > than gene_seq size.
                if (position_insertion > normalized_gene_seq.size()) {
                    std::string error_message = "";
                    error_message += "?Error: Insertion mutation index is greater than the length of the gene sequence.";
                    return error_message;
                }
                
                // Determine insertion sequence. Insertion sequence is substring that follows "ins" in mutation.
                std::string insertion_sequence = mutation.substr(mutation.find("ins") + 3);
                
                /* Generate lowercase insertion mutated gene sequence. */
                mutated_gene_seq = normalized_gene_seq.substr(0, position_insertion+1);
                mutated_gene_seq += to_lower_str(insertion_sequence);
                mutated_gene_seq += normalized_gene_seq.substr(position_insertion+1);
                
                // Insertion mutation successful.
                return mutated_gene_seq;
                
            } else {
                // Unknown/unsupported mutation
                std::string error_message = "?Error: Unknown mutation.";
                return error_message;
            }
            
            
        /* Unrecognized/Unsupported mutation strings. Error messages start with "?" to indicate error.
           int main() which utilizes nucleotide_utisl should catch these errors using "?" as keyword and should
            cut them out using .substr(1).
        */
        } else if (protein) {
            return "?Error: Protein reference sequence mutations are not yet supported.";
        } else if (rna) {
            return "?Error: RNA reference sequence mutations are not yet supported. Try converting to DNA mutation notation.";
        } else {
            return "?Error: Mutation was not recognized. Reference sequence type was not recognized. Examples of valid reference sequence types: c., g., m.";
        }
        
    }

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
    std::string extract_old_sequence(const std::string& gene_seq, const std::string& mutation){
        // Determine reference sequence.
        bool coding = mutation.find('c') != std::string::npos;
        bool genomic = mutation.find('g') != std::string::npos;
        bool mitochondrial = mutation.find('m') != std::string::npos;
        bool protein = mutation.find('p') != std::string::npos;
        bool rna = mutation.find('r') != std::string::npos;
        
        if (coding || genomic || mitochondrial) {   // All are DNA reference sequence mutations.
            
            // Determine mutation type.
            bool substitution = mutation.find('>') != std::string::npos;
            bool deletion = mutation.find("del") != std::string::npos;
            bool duplication = mutation.find("dup") != std::string::npos;
            bool insertion = mutation.find("ins") != std::string::npos;
            
            /* Determine given mutation sequence. */
            if (substitution) {
                /* Format of substitution mutation:
                    "c.123A>G"
                 */
                
                /* Determine old mucleotide
                    old_nucleotide will be to the left of '>'
                 */
                // Index of '>' within string mutation
                size_t index_greaterthan = mutation.find('>');
                size_t index_old_nucleotide = index_greaterthan - 1;
                char old_nucleotide = mutation[index_old_nucleotide];
                
                // Return old sequence
                return std::string (1, old_nucleotide);
                
            } else if (deletion) {
                /* Formats of deletion mutations:
                    "c.546delT"  --> single base at position 546 deleted.
                    "c.586_588del"  --> bases at positions 586, 587, 588 deleted.
                 */
                
                // Determine if single base deletion or multiple.
                bool multiple_deletions = mutation.find('_') != std::string::npos;
                if (multiple_deletions) {   // Multiple base deletion
                    /* Deal with indexing. */
                    // Generate substring of mutation to include only FIRST number.
                    std::string mutation_substr_1 = mutation.substr(0, mutation.find('_'));
                    // Generate substring of mutation to include only SECOND number.
                    std::string mutation_substr_2 = mutation.substr(mutation.find('_'));
                    // Find index of first nucleotide to be deleted.
                    size_t position_first = extract_num_from_str(mutation_substr_1) - 1; //0-BASED
                    // Find index of last nucleotide to be deleted.
                    size_t position_last = extract_num_from_str(mutation_substr_2) - 1; //0-BASED
                    // Calculate number of nucleotides to be deleted.
                    size_t num_deletions = position_last - position_first + 1;
                    
                    // Make sure given index is not > than gene_seq size.
                    if (position_last > gene_seq.size()) {
                        std::string error_message = "";
                        error_message += "?Error: Deletion mutation index is greater than the length of the gene sequence.";
                        return error_message;
                    }
                    
                    // Find old sequence (sequence to be deleted).
                    std::string old_sequence = gene_seq.substr(position_first, num_deletions);
                    
                    // Success.
                    return old_sequence;
                    
                } else {    // Single base deletion.
                    // Find deletion index.
                    // Subtract 1 to convert to 0-based indexing.
                    size_t position_mutation = extract_num_from_str(mutation) - 1; // 0-BASED
                    
                    // Make sure given index is not > than gene_seq size.
                    if (position_mutation > gene_seq.size()) {
                        std::string error_message = "";
                        error_message += "?Error: Deletion mutation index is greater than the length of the gene sequence.";
                        return error_message;
                    }
                    
                    // Find nucleotide to-be-deleted.
                    char old_nucleotide = gene_seq[position_mutation];
                    
                    // Success.
                    return std::string (1, old_nucleotide);
                }
                
            } else if (duplication) {
                /* Formats of duplication mutation:
                    "c.546dupT"  --> duplicate the nucleotide at position 546.
                    "c.586_588dup"  --> duplicate the nucleotide sequence 586, 587, 588.
                 */
                
                // Determine if single base or sequence duplication.
                bool sequence_duplication = mutation.find('_') != std::string::npos;
                
                /* Identify variables. */
                size_t position_first;
                size_t position_last;
                size_t length_duplication_sequence;
                
                if (!sequence_duplication) { // Single base duplication.
                    /* Deal with indexing for single base duplication.*/
                    // Find index of nucleotide in sequence-to-be-duplicated.
                    position_first = extract_num_from_str(mutation) - 1; // 0-BASED
                    position_last = position_first;
                    
                    // 1 nucleotide to be duplicated.
                    length_duplication_sequence = 1;
                    
                } else { // Sequence duplication.
                    /* Deal with indexing for sequence duplication.*/
                    // Generate substring of mutation to include only FIRST number.
                    std::string mutation_substr_1 = mutation.substr(0, mutation.find('_'));
                    // Generate substring of mutation to include only SECOND number.
                    std::string mutation_substr_2 = mutation.substr(mutation.find('_'));
                    // Find index of first nucleotide in sequence-to-be-duplicated.
                    position_first = extract_num_from_str(mutation_substr_1) - 1; // 0-BASED
                    // Find index of last nucleotide in sequence-to-be-duplicated.
                    position_last = extract_num_from_str(mutation_substr_2) - 1; // 0-BASED
                    // Calculate number of nucleotides to be duplicated.
                    length_duplication_sequence = position_last - position_first + 1;
                }
                
                
                // Generate sequence-to-be-duplicated.
                std::string sequence_to_be_duplicated = gene_seq.substr(position_first, length_duplication_sequence);
                                    
                // Return sequence-to-be-duplicated.
                return sequence_to_be_duplicated;
                
            } else if (insertion) {
                /* Format of insertion mutation:
                    "c.546_547insT"  --> insertion of T between 546 and 547
                    "c.186_187insGCGTGA"  --> insertion of GCGTGA between 186 and 187
                 */
                
                // Always return empty string.
                return "";
                
            } else {
                // Unknown/unsupported mutation
                std::string error_message = "?Error: Unknown mutation.";
                return error_message;
            }
            
            
        /* Unrecognized/Unsupported mutation strings. Error messages start with "?" to indicate error.
           int main() which utilizes nucleotide_utisl should catch these errors using "?" as keyword and should
            cut them out using .substr(1).
         */
        } else if (protein) {
            return "?Error: Protein reference sequence mutations are not yet supported.";
        } else if (rna) {
            return "?Error: RNA reference sequence mutations are not yet supported. Try converting to DNA mutation notation.";
        } else {
            return "?Error: Mutation was not recognized. Reference sequence type was not recognized. Examples of valid reference sequence types: c., g., m.";
        }
        
    }

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
    std::string extract_mutation_sequence(const std::string& gene_seq, const std::string& mutation){
        // Determine reference sequence.
        bool coding = mutation.find('c') != std::string::npos;
        bool genomic = mutation.find('g') != std::string::npos;
        bool mitochondrial = mutation.find('m') != std::string::npos;
        bool protein = mutation.find('p') != std::string::npos;
        bool rna = mutation.find('r') != std::string::npos;
        
        if (coding || genomic || mitochondrial) {   // All are DNA reference sequence mutations.
            
            // Determine mutation type.
            bool substitution = mutation.find('>') != std::string::npos;
            bool deletion = mutation.find("del") != std::string::npos;
            bool duplication = mutation.find("dup") != std::string::npos;
            bool insertion = mutation.find("ins") != std::string::npos;
            
            /* Determine given mutation sequence. */
            if (substitution) {
                /* Format of substitution mutation:
                    "c.123A>G"
                 */
                
                /* Determine mutated mucleotide
                    mutated_nucleotide will be to the right of '>'
                 */
                // Index of '>' within string mutation
                size_t index_greaterthan = mutation.find('>');
                size_t index_mutated_nucleotide = index_greaterthan + 1;
                char mutated_nucleotide = mutation[index_mutated_nucleotide];
                
                // Return mutated sequence
                return std::string (1, mutated_nucleotide);
                
            } else if (deletion) {
                // Deletion always returns empty string ""
                return "";
                
            } else if (duplication) {
                /* Formats of duplication mutation:
                    "c.546dupT"  --> duplicate the nucleotide at position 546.
                    "c.586_588dup"  --> duplicate the nucleotide sequence 586, 587, 588.
                 */
                
                // Determine if single base or sequence duplication.
                bool sequence_duplication = mutation.find('_') != std::string::npos;
                
                /* Identify variables. */
                size_t position_first;
                size_t position_last;
                size_t length_duplication_sequence;
                
                if (!sequence_duplication) { // Single base duplication.
                    /* Deal with indexing for single base duplication.*/
                    // Find index of nucleotide in sequence-to-be-duplicated.
                    position_first = extract_num_from_str(mutation) - 1; // 0-BASED
                    position_last = position_first;
                    
                    // 1 nucleotide to be duplicated.
                    length_duplication_sequence = 1;
                    
                } else { // Sequence duplication.
                    /* Deal with indexing for sequence duplication.*/
                    // Generate substring of mutation to include only FIRST number.
                    std::string mutation_substr_1 = mutation.substr(0, mutation.find('_'));
                    // Generate substring of mutation to include only SECOND number.
                    std::string mutation_substr_2 = mutation.substr(mutation.find('_'));
                    // Find index of first nucleotide in sequence-to-be-duplicated.
                    position_first = extract_num_from_str(mutation_substr_1) - 1; // 0-BASED
                    // Find index of last nucleotide in sequence-to-be-duplicated.
                    position_last = extract_num_from_str(mutation_substr_2) - 1; // 0-BASED
                    // Calculate number of nucleotides to be duplicated.
                    length_duplication_sequence = position_last - position_first + 1;
                }
                
                
                // Generate sequence-to-be-duplicated.
                std::string sequence_to_be_duplicated = gene_seq.substr(position_first, length_duplication_sequence);
                
                // Generate duplicated mutation sequence.
                std::string duplicated_seq = sequence_to_be_duplicated + sequence_to_be_duplicated;
                    
                // Return duplicated sequence.
                return duplicated_seq;
                
            } else if (insertion) {
                /* Format of insertion mutation:
                    "c.546_547insT"  --> insertion of T between 546 and 547
                    "c.186_187insGCGTGA"  --> insertion of GCGTGA between 186 and 187
                 */
                                
                // Determine insertion sequence. Insertion sequence is substring that follows "ins" in mutation.
                std::string insertion_sequence = mutation.substr(mutation.find("ins") + 3);
                
                // Return intertion sequence.
                return insertion_sequence;
                
            } else {
                // Unknown/unsupported mutation
                std::string error_message = "?Error: Unknown mutation.";
                return error_message;
            }
            
            
        /* Unrecognized/Unsupported mutation strings. Error messages start with "?" to indicate error.
           int main() which utilizes nucleotide_utisl should catch these errors using "?" as keyword and should
            cut them out using .substr(1).
         */
        } else if (protein) {
            return "?Error: Protein reference sequence mutations are not yet supported.";
        } else if (rna) {
            return "?Error: RNA reference sequence mutations are not yet supported. Try converting to DNA mutation notation.";
        } else {
            return "?Error: Mutation was not recognized. Reference sequence type was not recognized. Examples of valid reference sequence types: c., g., m.";
        }
        
    }

    /*
     * Determines and returns the 0-BASED index of the FIRST nucleotide edited by a given mutation.
     * param mutation is given in standard mutation notation (see function comment nucleotide::mutate()).
     * Ex. extract_mutation_position("c.586_591dup") returns 567 (0-based indexing subtracts 1 from 1-based)
     */
    size_t extract_mutation_position(const std::string& mutation){
        return extract_num_from_str(mutation) - 1;
    }

    /*
     * Estimates the "likelihood" [0.0, 1.0] that two nucleic acid sequences will form a stable
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
    double dimerize(std::string seq1, std::string seq2) {
        auto to_dna = [](std::string s) {
            for (char& c : s) {
                c = (char)std::toupper((unsigned char)c);
                if (c == 'U') c = 'T';
            }
            return s;
        };

        auto comp = [](char b) -> char {
            switch (b) {
                case 'A': return 'T';
                case 'T': return 'A';
                case 'C': return 'G';
                case 'G': return 'C';
                default:  return '?';
            }
        };

        auto revcomp = [&](const std::string& s) {
            std::string rc;
            rc.reserve(s.size());
            for (auto it = s.rbegin(); it != s.rend(); ++it) {
                char c = comp(*it);
                rc.push_back(c);
            }
            return rc;
        };

        auto is_wc_comp = [&](char a, char b) {
            // a on strand1, b on opposite strand (in same index alignment)
            if (a=='A' && b=='T') return true;
            if (a=='T' && b=='A') return true;
            if (a=='C' && b=='G') return true;
            if (a=='G' && b=='C') return true;
            return false;
        };

        auto is_gc_pair = [&](char a, char b) {
            return (a=='G' && b=='C') || (a=='C' && b=='G');
        };

        seq1 = to_dna(seq1);
        seq2 = to_dna(seq2);

        // Align seq1 vs reverse-complement of seq2 (the physical duplex orientation).
        const std::string rc2 = revcomp(seq2);

        // SantaLucia 1998 ΔG°37 nearest-neighbor stack parameters (kcal/mol)
        // Format: "XY/ZW" where XY is 5'->3' on seq1, ZW is 3'->5' on opposite strand.
        const std::unordered_map<std::string, double> dG = {
            {"AA/TT", -1.00}, {"AT/TA", -0.88}, {"TA/AT", -0.58}, {"CA/GT", -1.45}, {"GT/CA", -1.44},
            {"CT/GA", -1.28}, {"GA/CT", -1.30}, {"CG/GC", -2.17}, {"GC/CG", -2.24}, {"GG/CC", -1.84}
        };

        auto stack_dG37_for_segment = [&](int i1, int i2, int len) -> double {
            // i1 indexes seq1, i2 indexes rc2, and len is a contiguous perfectly complementary segment length.
            // Sum NN stacks across the segment (len-1 stacks).
            double sum = 0.0;
            for (int k = 0; k < len - 1; ++k) {
                char a1 = seq1[i1 + k];
                char a2 = seq1[i1 + k + 1];

                // Opposite strand is rc2 (5'->3'), but NN keys want opposite dinuc in 3'->5' direction.
                char b1_5to3 = rc2[i2 + k];
                char b2_5to3 = rc2[i2 + k + 1];
                // Reverse it to get 3'->5' dinucleotide:
                std::string opp_3to5;
                opp_3to5.push_back(b2_5to3);
                opp_3to5.push_back(b1_5to3);

                std::string key;
                key.push_back(a1);
                key.push_back(a2);
                key.push_back('/');
                key += opp_3to5;

                auto it = dG.find(key);
                if (it == dG.end()) {
                    // If something unexpected happens (non-ATGC), treat as neutral (break strength).
                    return 0.0;
                }
                sum += it->second;
            }
            return sum;
        };

        // Scan all relative offsets between seq1 and rc2 to find:
        // - strongest (most negative) ΔG segment
        // - longest contiguous complementarity anywhere
        // - 3'-end complementarity within last 8 nt for each oligo (PCR-primer-style risk)
        double best_dG = 0.0;          // more negative = stronger
        int best_len_any = 0;

        int best_3p_len = 0;
        int best_3p_gc_len = 0;

        const int n1 = (int)seq1.size();
        const int n2 = (int)rc2.size();

        // offset = i2 - i1. i1 in [0,n1), i2 in [0,n2)
        for (int offset = -n1; offset <= n2; ++offset) {
            int i1 = std::max(0, -offset);
            int i2 = std::max(0,  offset);
            int L  = std::min(n1 - i1, n2 - i2);
            if (L <= 0) continue;

            int run = 0;
            int run_start1 = i1;
            int run_start2 = i2;

            for (int t = 0; t < L; ++t) {
                char a = seq1[i1 + t];
                char b = rc2[i2 + t];

                bool ok = (a=='A'||a=='T'||a=='C'||a=='G') && (b=='A'||b=='T'||b=='C'||b=='G') && is_wc_comp(a,b);

                if (ok) {
                    if (run == 0) {
                        run_start1 = i1 + t;
                        run_start2 = i2 + t;
                    }
                    run++;

                    // track GC-run length within this run as well
                } else {
                    if (run >= 2) {
                        best_len_any = std::max(best_len_any, run);
                        double seg_dG = stack_dG37_for_segment(run_start1, run_start2, run);
                        if (seg_dG < best_dG) best_dG = seg_dG;

                        // 3' end emphasis: does this segment touch within last 8 nt of BOTH oligos?
                        // 3' end of seq1 is index n1-1, 3' end of rc2 corresponds to 3' end of seq2.
                        int end1 = run_start1 + run - 1;
                        int end2 = run_start2 + run - 1;
                        bool near3p_1 = (end1 >= n1 - 8);
                        bool near3p_2 = (end2 >= n2 - 8);
                        if (near3p_1 && near3p_2) {
                            best_3p_len = std::max(best_3p_len, run);

                            // compute best contiguous GC-only within this segment
                            int gc_run = 0, best_gc = 0;
                            for (int k = 0; k < run; ++k) {
                                char aa = seq1[run_start1 + k];
                                char bb = rc2[run_start2 + k];
                                if (is_gc_pair(aa, bb)) { gc_run++; best_gc = std::max(best_gc, gc_run); }
                                else gc_run = 0;
                            }
                            best_3p_gc_len = std::max(best_3p_gc_len, best_gc);
                        }
                    }
                    run = 0;
                }
            }

            // flush run at end
            if (run >= 2) {
                best_len_any = std::max(best_len_any, run);
                double seg_dG = stack_dG37_for_segment(run_start1, run_start2, run);
                if (seg_dG < best_dG) best_dG = seg_dG;

                int end1 = run_start1 + run - 1;
                int end2 = run_start2 + run - 1;
                bool near3p_1 = (end1 >= n1 - 8);
                bool near3p_2 = (end2 >= n2 - 8);
                if (near3p_1 && near3p_2) {
                    best_3p_len = std::max(best_3p_len, run);

                    int gc_run = 0, best_gc = 0;
                    for (int k = 0; k < run; ++k) {
                        char aa = seq1[run_start1 + k];
                        char bb = rc2[run_start2 + k];
                        if (is_gc_pair(aa, bb)) { gc_run++; best_gc = std::max(best_gc, gc_run); }
                        else gc_run = 0;
                    }
                    best_3p_gc_len = std::max(best_3p_gc_len, best_gc);
                }
            }
        }

        // Convert features -> probability-ish score in [0,1]
        // IDT guidance often flags dimers more concerning when ΔG is more negative than about -9 kcal/mol
        auto logistic = [](double x) { return 1.0 / (1.0 + std::exp(-x)); };

        // If best_dG is 0, there was no ≥2-bp contiguous WC segment; risk should be very low.
        // Map ΔG: more negative => higher risk.
        // Center around -9, slope ~ 1.5 kcal/mol.
        double dg_score = logistic(((-9.0) - best_dG) / 1.5);  // e.g., best_dG=-12 => positive => high

        // Map 3' complementarity: 4+ contiguous at 3' ends = danger.
        double threep_score = 0.0;
        if (best_3p_len >= 4) threep_score = 0.8;
        else if (best_3p_len == 3) threep_score = 0.5;
        else if (best_3p_len == 2) threep_score = 0.2;

        // Extra boost for GC-rich 3' pairing (stronger binding).
        double threep_gc_boost = 0.0;
        if (best_3p_gc_len >= 3) threep_gc_boost = 0.2;
        else if (best_3p_gc_len == 2) threep_gc_boost = 0.1;

        // Long perfect complementarity anywhere also matters (non-3' dimers/sequestration).
        double any_score = 0.0;
        if (best_len_any >= 8) any_score = 0.6;
        else if (best_len_any >= 6) any_score = 0.4;
        else if (best_len_any >= 4) any_score = 0.2;

        // Combine (bounded to [0,1])
        double p = 0.55 * dg_score + 0.30 * (threep_score + threep_gc_boost) + 0.15 * any_score;
        p = std::clamp(p, 0.0, 1.0);
        return p;
    }

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
    double hairpin(std::string seq) {
        /* Tunable constants */
        constexpr std::size_t MIN_LOOP = 3;
        constexpr std::size_t MAX_LOOP = 30;
        constexpr std::size_t MIN_STEM_BP = 3;
        constexpr std::size_t THREE_PRIME_WINDOW = 8;

        auto to_dna = [](std::string s) {
            for (char& c : s) {
                c = (char)std::toupper((unsigned char)c);
                if (c == 'U') c = 'T';
            }
            return s;
        };

        auto is_ATGC = [](char c) {
            return c=='A' || c=='T' || c=='C' || c=='G';
        };

        auto comp = [](char b) -> char {
            switch (b) {
                case 'A': return 'T';
                case 'T': return 'A';
                case 'C': return 'G';
                case 'G': return 'C';
                default:  return '?';
            }
        };

        auto is_wc_comp = [&](char a, char b) {
            return is_ATGC(a) && is_ATGC(b) && (a == comp(b));
        };

        seq = to_dna(seq);
        const std::size_t N = seq.size();
        if (N < (2 * MIN_STEM_BP + MIN_LOOP)) return 0.0;

        // SantaLucia 1998 ΔG°37 nearest-neighbor stack parameters (kcal/mol)
        // Format: "XY/ZW" where XY is 5'->3' on one strand, ZW is 3'->5' on the opposite.
        const std::unordered_map<std::string, double> dG = {
            {"AA/TT", -1.00}, {"AT/TA", -0.88}, {"TA/AT", -0.58}, {"CA/GT", -1.45}, {"GT/CA", -1.44},
            {"CT/GA", -1.28}, {"GA/CT", -1.30}, {"CG/GC", -2.17}, {"GC/CG", -2.24}, {"GG/CC", -1.84}
        };

        // Score a perfectly complementary stem of length 'len' pairing:
        //   left segment:  seq[l0 ... l0+len-1]   (5'->3')
        //   right segment: seq[r0 ... r0+len-1]   (5'->3', downstream)
        // Base pairing: seq[l0+k] pairs with seq[r0+(len-1-k)] (antiparallel in hairpin stem).
        auto stem_stack_dG37 = [&](std::size_t l0, std::size_t r0, std::size_t len) -> double {
            if (len < 2) return 0.0;
            double sum = 0.0;

            // Each stack corresponds to adjacent base pairs in the stem.
            for (std::size_t k = 0; k + 1 < len; ++k) {
                char L1 = seq[l0 + k];
                char L2 = seq[l0 + k + 1];

                // Paired bases on the right (5'->3' indexing), but traversed in 3'->5' order
                char R1 = seq[r0 + (len - 1 - k)];
                char R2 = seq[r0 + (len - 2 - k)];

                if (!is_ATGC(L1) || !is_ATGC(L2) || !is_ATGC(R1) || !is_ATGC(R2)) return 0.0;
                if (!is_wc_comp(L1, R1) || !is_wc_comp(L2, R2)) return 0.0;

                std::string key;
                key.push_back(L1);
                key.push_back(L2);
                key.push_back('/');
                key.push_back(R1);
                key.push_back(R2);

                auto it = dG.find(key);
                if (it == dG.end()) return 0.0;
                sum += it->second;
            }
            return sum;
        };

        // Track strongest hairpin candidate
        double best_dG = 0.0;            // most negative = strongest
        std::size_t best_stem = 0;       // longest contiguous stem found
        std::size_t best_3p_stem = 0;    // longest stem that involves the 3' end window
        std::size_t best_3p_gc_run = 0;  // GC run within best 3' stem (simple proxy)

        // Enumerate hairpin geometries:
        // For a chosen loop length and left-stem start l0, a stem of length L implies right-stem start:
        //   r0 = l0 + L + loop
        for (std::size_t loop = MIN_LOOP; loop <= MAX_LOOP; ++loop) {
            if (2 * MIN_STEM_BP + loop > N) break;

            for (std::size_t l0 = 0; l0 + 2 * MIN_STEM_BP + loop <= N; ++l0) {
                const std::size_t max_len = (N - l0 - loop) / 2;
                if (max_len < MIN_STEM_BP) continue;

                // Find the longest valid stem for this (l0, loop) by checking candidate L from max down.
                for (std::size_t L = max_len; L >= MIN_STEM_BP; --L) {
                    const std::size_t r0 = l0 + L + loop;
                    if (r0 + L > N) {
                        if (L == MIN_STEM_BP) break;
                        continue;
                    }

                    bool ok = true;
                    for (std::size_t k = 0; k < L; ++k) {
                        const char left  = seq[l0 + k];
                        const char right = seq[r0 + (L - 1 - k)];
                        if (!is_wc_comp(left, right)) { ok = false; break; }
                    }

                    if (!ok) {
                        if (L == MIN_STEM_BP) break;
                        continue;
                    }

                    // Longest stem found for this (l0, loop).
                    best_stem = std::max(best_stem, L);

                    const double seg_dG = stem_stack_dG37(l0, r0, L);
                    if (seg_dG < best_dG) best_dG = seg_dG;

                    // 3' involvement: if either side overlaps last THREE_PRIME_WINDOW bases
                    const std::size_t threep_start = (N > THREE_PRIME_WINDOW) ? (N - THREE_PRIME_WINDOW) : 0;
                    const bool touches_3p =
                        (l0 + L - 1 >= threep_start) || (r0 + L - 1 >= threep_start);

                    if (touches_3p) {
                        best_3p_stem = std::max(best_3p_stem, L);

                        std::size_t gc_run = 0, best_gc = 0;
                        for (std::size_t k = 0; k < L; ++k) {
                            const char c = seq[l0 + k];
                            if (c == 'G' || c == 'C') { ++gc_run; best_gc = std::max(best_gc, gc_run); }
                            else gc_run = 0;
                        }
                        best_3p_gc_run = std::max(best_3p_gc_run, best_gc);
                    }

                    break; // stop at the longest L for this (l0, loop)
                }
            }
        }

        // Feature -> likelihood mapping
        auto logistic = [](double x) { return 1.0 / (1.0 + std::exp(-x)); };

        // Hairpins can be relevant at less negative ΔG than primer dimers (intramolecular is entropically favored).
        const double dg_score = logistic(((-6.5) - best_dG) / 1.2);

        double len_score = 0.0;
        if      (best_stem >= 10) len_score = 0.75;
        else if (best_stem >=  8) len_score = 0.55;
        else if (best_stem >=  6) len_score = 0.35;
        else if (best_stem >=  4) len_score = 0.20;
        else if (best_stem >=  3) len_score = 0.10;

        double threep_score = 0.0;
        if      (best_3p_stem >= 6) threep_score = 0.30;
        else if (best_3p_stem >= 4) threep_score = 0.20;
        else if (best_3p_stem >= 3) threep_score = 0.10;

        double threep_gc_boost = 0.0;
        if      (best_3p_gc_run >= 4) threep_gc_boost = 0.10;
        else if (best_3p_gc_run >= 3) threep_gc_boost = 0.06;
        else if (best_3p_gc_run >= 2) threep_gc_boost = 0.03;

        const double p = 0.55 * dg_score + 0.30 * len_score + 0.15 * (threep_score + threep_gc_boost);
        return std::clamp(p, 0.0, 1.0);
    }


} // end namespace nucleotide

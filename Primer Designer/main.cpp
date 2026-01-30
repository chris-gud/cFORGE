/*
 * primer-designer
 *
 * A C++ program for designing PCR primers for gene cloning, site-directed
 * mutagenesis, and gene insertion into plasmid vectors. The program operates
 * locally, parses standard FASTA sequence files and JSON-formatted restriction
 * enzyme libraries, and generates validated primer sets along with optional
 * plasmid map visualizations.
 *
 * 2026 Chris Gudmundsen (https://github.com/chris-gud)
 * Released under the MIT License. See LICENSE for details.
 */

#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <stdexcept>
#include <utility>

#include "Primer.hpp"
#include "RestrictionEnzyme.hpp"
#include "SequenceRecord.hpp"
#include "nucleotide_utils.hpp"
#include "json.hpp"


// Filepath for restriction enzyme JSON file.
// User can edit if desires specific restriction enzymes.
// Use restriction_enzyme_json repo on my Github to create a json file with the correct attributes.
const std::string RE_sites_file_path = "data/RE_sites1.json";

/* All may not be used depending on user instruction. */
std::string vector_filename, gene_filename;
std::string vector_seq, vector_header, gene_seq, gene_header;
std::string relative_path = "data/";


/* --------------------------------------------------------------------------------------------- */
/* ------------------------------------------- HEADERS ----------------------------------------- */
// Function comments are located above function definitions.
bool read_fasta(const std::string& filename, std::string& sequence, std::string& header);

std::vector<RestrictionEnzyme> read_RE_json(const std::string& filename); // Should put this in RE file

size_t determine_cut_index(const size_t& find_index, const size_t& local_RE_cut_index);

bool found_again(const size_t& index_found, const std::string& recognition_seq, const std::string& main_seq);

std::string find_available_filename(const std::string& base_name, const std::string& extension);

std::string extract_vector_name(const std::string& fasta_header);

void write_linear_cut_map(
    const std::string vector_name,
    const std::string &vector_seq,
    const std::vector<RestrictionEnzyme> &enzymes
);

void write_svg_vector_map(
    const std::string &vector_name,
    const std::string &vector_seq,
    const std::vector<RestrictionEnzyme> &enzymes
);

std::string query_user (
    const std::string& query,
    const std::vector<std::string>& valid_inputs,
    const std::string& invalid_input_message
);

std::string repeat(std::string str, size_t num_times);

std::string generate_tag_seq(std::string& tag);

std::string generate_protease_seq(std::string& protease_id);


/* ----------------------------------------- END HEADERS --------------------------------------- */
/* --------------------------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------------------------- */
/* --------------------------------------------- MAIN ------------------------------------------ */
int main() {
    /* Variables for use any time a user is queried throughout the program. */
    std::string output;
    std::vector<std::string> valid_inputs;
    std::string error_message;
    std::string user_input;
    
    /* Ask which task user wants to perform. */
    output = "Functions:"
             "\n  A.) Generate flanking primers for gene cloning."
             "\n  B.) Generate mutagenic primers for PCR mutagenesis."
             "\n  C.) Generate primers for insertion of a gene into a vector."
             "\n  D.) Check for dimerization between two primers."
             "\n  E.) Generate the reverse complement of a sequence."
             "\n\nEnter desired function";
    valid_inputs = {"A", "a", "flanking", "B", "b", "mutagenic", "C", "c", "insertion", "D", "d", "dimerization", "E", "e", "revcomp"};
    error_message = "That option was not recognized. Try again.";
    user_input = query_user(output, valid_inputs, error_message);
    
    /* --------------------------------------------------------------------------------------------- */
    /*--------------------------------- FLANKING PRIMERS ------------------------------------------- */

    if (user_input == "A" || user_input == "a" || user_input == "flanking") {
        /* A.) GENERATE FLANKING PRIMERS FOR GENE CLONING. */
        std::cout << "\n- GENERATING FLANKING PRIMERS FOR GENE CLONING -";
        
        // Read in gene FASTA only
        std::cout << "\nEnter gene FASTA filename: ";
        std::getline(std::cin, gene_filename);
        
        if (!read_fasta(relative_path + gene_filename, gene_seq, gene_header)) {
            std::cerr << "Error reading gene file.\n";
            return 1;
        }
        
        if (gene_seq.size() < 36) {
            std::cerr << "Gene sequence " << gene_seq << " is too small.\n";
            return 1;
        }
        
        // Generate primers
        std::vector<Primer> forward_primers; // Will be ordered based on length
        std::vector<Primer> reverse_primers; // "
        
        // Generate list of FORWARD primers (forward_primers)
        for (size_t i = 18; i <= 25; i++) {
            // Forward primers bind to the reverse strand
            std::string primer_seq = gene_seq.substr(0, i);
            // Create a potential primer
            Primer x = Primer(primer_seq);
            // Determine if x satisfies all internal rules
            if (x.is_good()) {
                // Ensure x doesn't bind off-target
                if (!x.binds_off_target(gene_seq)) {
                    // All internal rules met + Primer does not bind off-target; x added to list of viable FORWARD primers
                    forward_primers.push_back(x);
                }
            }
        }
        
        // Generate list of REVERSE primers (reverse_primers)
        for (size_t i = 18; i <= 25; i++) {
            // Note: Reverse primers bind to the forward strand.
            std::string forward_end = gene_seq.substr(gene_seq.size() - i);
            // Reverse complement of the end of the forward seqeunce of the gene is the reverse primer.
            std::string reverse_primer_seq = nucleotide::revcomp_of(forward_end);
            // Create a potential primer
            Primer x = Primer(reverse_primer_seq);
            // Determine if x satisfies all internal rules
            if (x.is_good()) {
                // Ensure x doesn't bind off-target
                if (!x.binds_off_target(gene_seq)) {
                    // All internal rules met + Primer does not bind off-target; x added to list of viable REVERSE primers
                    reverse_primers.push_back(x);
                }
            }
        }
        
        
        
        // Generate .txt file containing the approved forward and reverse primers
        // Find an available filename
        std::string primer_filename = find_available_filename("output/cloning_primers", ".txt");
        // Create and open a text file for writing in the output folder
        std::ofstream outputFile(primer_filename);
        // Check if the file was successfully opened
        if (outputFile.is_open()) {
            // Write data to file
            outputFile << " - VIABLE PRIMERS FOR GENE CLONING -\n\n"
                       << "Gene source file: " << gene_filename << "\n";
            
            // Print FORWARD primers to terminal and write to output file
            std::string for_header = "\n"
                                     "|------------------------------------------|\n"
                                     "| Viable FORWARD Primers for gene cloning: |\n"
                                     "|------------------------------------------|\n";
            outputFile << for_header;
            for (Primer p : forward_primers) {
                std::string out = p.to_string();
                //std::cout << out;
                outputFile << out;
            }
            
            // Print REVERSE primers to terminal and write to output file
            std::string rev_header = "\n\n"
                                     "|------------------------------------------|\n"
                                     "| Viable REVERSE Primers for gene cloning: |\n"
                                     "|------------------------------------------|\n";
            outputFile << rev_header;
            for (Primer r : reverse_primers) {
                std::string out = r.to_string();
                //std::cout << out;
                outputFile << out;
            }
            
            // Close  file
            outputFile.close();

            std::cout << "\nSuccess: A list of viable primers has been output to " << primer_filename.substr(7)
                      << "\nSelect your two primers and ensure they do not dimerize by using function (D) of this program.\n\n";
            
            // Success
            return 0;
        } else {
            std::cerr << "\nError creating primer output file: " << primer_filename << "\n";
            return 1;
        }
    } 
    /*-------------------------------- END FLANKING PRIMERS ---------------------------------------- */
    /* --------------------------------------------------------------------------------------------- */

    /* --------------------------------------------------------------------------------------------- */
    
    /* --------------------------------------------------------------------------------------------- */
    /*--------------------------------- MUTAGENIC PRIMERS ------------------------------------------ */
    
    else if (user_input == "B" || user_input == "b" || user_input == "mutagenic") {
        /* B.) GENERATE MUTAGENIC PRIMERS FOR PCR MUTAGENESIS. */
        std::cout << "\n- GENERATING MUTAGENIC PRIMERS FOR PCR MUTAGENESIS -";
        
        std::cout << "\nEnter gene FASTA filename: ";
        std::getline(std::cin, gene_filename);
                
        if (!read_fasta(relative_path + gene_filename, gene_seq, gene_header)) {
            std::cerr << "Error reading gene file.\n";
            return 1;
        }
        
        /* Loop until user input mutation is read. */
        std::string user_mutation;
        std::string old_sequence;
        std::string mutated_gene_seq;
        std::string mutation_sequence;
        size_t mutation_start_index = 0;
        size_t mutation_end_index = 0;

        bool is_good_input = false;
        while (!is_good_input) {
            // Query user which position to mutate.
            std::cout << "\nEnter desired mutation in standard notation: ";
            std::getline(std::cin, user_mutation);
            
            // If escape string ("q"), quit program.
            if (user_mutation == "q") {
                std::cout << "\nQuitting program.\n";
                break;
            }
            
            // Generate mutated gene sequence.
            mutated_gene_seq = nucleotide::mutate(gene_seq, user_mutation);
                
            /* Check for error code (string returned by nucleotide::mutate() begins with '?').
                If error code, continue;
             */
            if (mutated_gene_seq[0] == '?') {
                std::cout << "\n" << mutated_gene_seq.substr(1) << "\nPlease try again. Enter 'q' to exit.\n";
                continue;
            }
            
            /* Assign values to variables. */
            old_sequence = nucleotide::extract_old_sequence(gene_seq, user_mutation);
            mutation_sequence = nucleotide::extract_mutation_sequence(gene_seq, user_mutation);
            mutation_start_index = nucleotide::extract_mutation_position(user_mutation);
            mutation_end_index = mutation_start_index + mutation_sequence.size(); // Non-inclusive.
            
            /* Check if mutation is within 13 nucleotides of start or end of gene.
               If mutation is within 13 nucleotides, PCR mutagenesis is likely not viable for this type of mutation.
                Continue next iteration of loop.
               If mutation is NOT within 13 nucleotides, input is good.
                good_input = true. Finish loop then break.
             */
            if (mutation_start_index <= 13 || gene_seq.size() - mutation_end_index <= 13) {
                std::cout << "\nThe given mutation is within 13 nucleotides of the start or end of your gene; site directed PCR by mutagenesis is not a viable method to confer this mutation. Please try again. Enter 'q' to exit.\n";
                continue;
            } else {
                is_good_input = true;
                break;
            }
        }
        
        std::cout << "\nMutation Position = " << mutation_start_index + 1;
        std::cout << "\nOld Sequence = \"" << old_sequence << "\"";
        std::cout << "\nNew Mutated Sequence = \"" << mutation_sequence << "\"\n";

        
        
        /* Make pseudo-primer Primer objects.
           Pseudo-primers are sequences that consist of only the flanking homology that anneals to the gene DNA; the mutated nucleotide(s) is/are not included.
            Ex.
                Mutagenic reverse primer sequence:         AAAAAAAAAAtAAAAAAAAA
                Gene sequence:                       ...GGGTTTTTTTTTTGTTTTTTTTTCCC...
                --> Reverse Pseudo-primer sequence:        AAAAAAAAAA-AAAAAAAAA
                
                The mutation is the G->T substitution. The mutated T does not bind to the original G; therefore it does not contribute to the melting temp, GC-content, max length, etc. For that reason, we create a pseudo-primer Primer object consisting of only the annealing homology to utilize functions within the Primer object, then we generate the real primer sequences for the user afterwards.
         
                NOTE: there is a systematic error involved in this workflow. Some of the functionality involved in the Primer.is_good() method is not transferable between pseudo and real primers. For example, the Primer.is_good() function of the Primer object will asses whether the pseudo-primer homodimerizes or forms a hairpin, and, if it does either, will not validate it. However, if the site of the would-be mutation contributes to homodimerization/hairpins (i.e. the absence of mutated sequence leads to homodimerization/hairpins), the is_good() method will return a spurious "false" flag. While this is obviously a mistake, in my opinion this is so unlikely to happen that it is not worth overwriting several methods or making a new Primer subclass to avoid this.
         */
        // Vectors of pseudo-primer Primer objects.
        std::vector<Primer> forward_pseudo_primers; // Will be ordered based on length
        std::vector<Primer> reverse_pseudo_primers; // "
        
        // Vectors of real primer string sequences.
        std::vector<std::string> forward_real_primers_seq; // Will be ordered based on length
        std::vector<std::string> reverse_real_primers_seq; // "
        
        // Generate list of viable primers.
        for (size_t i = 18; i <= 45; i++) {
            size_t half_i = i / 2;  // Will be rounded down.
            //size_t remainder = i % 2;
            
            // Calculate primer start index within the gene sequence.
            size_t primer_start_index = mutation_start_index - half_i;
            // Calculate length of primer + mutated sequence.
            size_t length_real_primer = i + mutation_sequence.size();
            
            /* Make sure proposed primer is not outside of scope of gene_seq. */
            // Note: primer_start_index is unsigned AND the mutation location is guaranteed to not be greater than the
            //  gene sequence length (by nucloetide_utils); the >= operator in this case ONLY checks for wraparound.
            if (primer_start_index >= gene_seq.size()) {
                primer_start_index = 0;
            }
            // Flagged if end of primer is greater than the gene size.
            if (primer_start_index + length_real_primer >= gene_seq.size()) {
                /* Offset the primer so its tail is right at the end of the gene sequence. */
                size_t offset = primer_start_index + length_real_primer - gene_seq.size();
                primer_start_index -= offset;
            }
            
            // Extract forward real primer string.
            std::string forward_real_primer_seq = mutated_gene_seq.substr(primer_start_index, length_real_primer);
            
            /* Determine lengths of each binding homology. */
            size_t length_binding_homology1 = mutation_start_index - primer_start_index;
            size_t length_binding_homology2 = i - length_binding_homology1;

            /* Generate pseudo-primer string. */
            std::string forward_pseudo_primer_seq = "";
            forward_pseudo_primer_seq += gene_seq.substr(primer_start_index, length_binding_homology1);
            forward_pseudo_primer_seq += gene_seq.substr(mutation_end_index, length_binding_homology2);
            
            // Create Primer object using pseudo-primer sequence.
            Primer forward_psuedo_primer = Primer(forward_pseudo_primer_seq);

            // Determine if forward_psuedo_primer satisfies all internal rules.
            if (forward_psuedo_primer.is_good()) {
                // Ensure forward_psuedo_primer doesn't bind off-target.
                if (!forward_psuedo_primer.binds_off_target(gene_seq)) {
                    // All internal rules met + Primer does not bind off-target; x added to list of viable FORWARD primers.
                    forward_pseudo_primers.push_back(forward_psuedo_primer); // Vector of Primers.
                    forward_real_primers_seq.push_back(forward_real_primer_seq); // Vector of std::strings.
                }
            }
            
            
            /* Generate reverse primer. Reverse mutagenic primers are the reverse complement of the forward mutagenic
                primers.
             */
            // Generate reverse primer string.
            std::string reverse_real_primer_seq = nucleotide::revcomp_of(forward_real_primer_seq);
            std::string reverse_pseudo_primer_seq = nucleotide::revcomp_of(forward_pseudo_primer_seq);
            
            // Generate reverse pseudo-primer Primer object.
            Primer reverse_pseudo_primer = Primer(reverse_pseudo_primer_seq);

            // Determine if reverse_pseudo_primer satisfies all internal rules.
            if (reverse_pseudo_primer.is_good()) {
                // Ensure reverse_pseudo_primer doesn't bind off-target.
                if (!reverse_pseudo_primer.binds_off_target(gene_seq)) {
                    // All internal rules met + Primer does not bind off-target; x added to list of viable REVERSE primers.
                    reverse_pseudo_primers.push_back(reverse_pseudo_primer); // Vector of Primers.
                    reverse_real_primers_seq.push_back(reverse_real_primer_seq); // Vector of std::strings.
                }
            }
        } // Lists of viable forward and reverse Primers have been generated.
        
        
        // Generate .txt file containing the approved forward and reverse primers for mutagenesis.
        // Find an available filename.
        std::string primer_filename = find_available_filename("output/mutagenic_primers", ".txt");
        // Create and open a text file for writing in the output folder.
        std::ofstream outputFile(primer_filename);
        // Check if the file was successfully opened.
        if (outputFile.is_open()) {
            // Write data to file.
            outputFile << " - VIABLE PRIMERS FOR SITE DIRECTED MUTAGENESIS BY PCR -\n\n"
                       << "Gene source file: " << gene_filename << "\n"
                       << "Mutation: " << user_mutation << "\n";
            
            // FORWARD primers write to output file.
            std::string forward_primer_title = "\n"
                                     "|---------------------------------------------|\n"
                                     "| Viable FORWARD Primers for PCR mutagenesis: |\n"
                                     "|---------------------------------------------|\n";
            outputFile << forward_primer_title;
            for (size_t i = 0; i < forward_pseudo_primers.size(); i++) {
                
                Primer* p = &forward_pseudo_primers[i];
                
                std::string string_seq = "5'- " + forward_real_primers_seq[i] + " -3'";
                std::string string_length = std::to_string(forward_real_primers_seq[i].size()) + " n.t.";
                std::string string_wallace_melting_temp = std::to_string(p->wallace_melting_temp) + "°C";
                std::string string_sl_melting_temp = std::to_string(p->santalucia_melting_temp).substr(0, 4) + "°C";
                std::string string_GC_content = std::to_string(p->GC_content*100.).substr(0, 4) + "%";
                std::string string_GC_clamp = p->GC_clamp ? "YES" : "NO";
                std::string string_homo_dimerizes = p->homo_dimerizes ? "YES" : "NO";
                std::string string_max_length_homopolymer_run = std::to_string(p->max_length_homopolymer_run) + " n.t.";
                std::string string_hairpin = p->hairpin ? "YES" : "NO";
               
                std::string out = "";
                out += "----------------------------------------------------------------------\n";
                out += "              Primer Sequence: " + string_seq + "\n";
                out += "                       Length: " + string_length + "\n";
                out += "         Wallace Melting Temp: " + string_wallace_melting_temp + "\n";
                out += "      SantaLucia Melting Temp: " + string_sl_melting_temp + "\n";
                out += " GC Content of Binding Region: " + string_GC_content + + "\n";
                out += "                     GC Clamp: " + string_GC_clamp + "\n";
                out += "             Homodimerization: " + string_homo_dimerizes + "\n";
                out += "       Homopolymer Run Length: " + string_max_length_homopolymer_run + "\n";
                out += "             Forms Hairpin(s): " + string_hairpin + "\n";
                out += "----------------------------------------------------------------------\n";
                
                //std::cout << out;
                outputFile << out;
            }
            
            // REVERSE primers write to output file
            std::string reverse_primer_title = "\n\n"
                                     "|---------------------------------------------|\n"
                                     "| Viable REVERSE Primers for PCR mutagenesis: |\n"
                                     "|---------------------------------------------|\n";
            outputFile << reverse_primer_title;
            for (size_t i = 0; i < reverse_pseudo_primers.size(); i++) {
                
                Primer* r = &reverse_pseudo_primers[i];
                
                std::string string_seq = "5'- " + reverse_real_primers_seq[i] + " -3'";
                std::string string_length = std::to_string(reverse_real_primers_seq[i].size()) + " n.t.";
                std::string string_wallace_melting_temp = std::to_string(r->wallace_melting_temp) + "°C";
                std::string string_sl_melting_temp = std::to_string(r->santalucia_melting_temp).substr(0, 4) + "°C";
                std::string string_GC_content = std::to_string(r->GC_content*100.).substr(0, 4) + "%";
                std::string string_GC_clamp = r->GC_clamp ? "YES" : "NO";
                std::string string_homo_dimerizes = r->homo_dimerizes ? "YES" : "NO";
                std::string string_max_length_homopolymer_run = std::to_string(r->max_length_homopolymer_run) + " n.t.";
                std::string string_hairpin = r->hairpin ? "YES" : "NO";
               
                std::string out = "";
                out += "----------------------------------------------------------------------\n";
                out += "              Primer Sequence: " + string_seq + "\n";
                out += "                       Length: " + string_length + "\n";
                out += "         Wallace Melting Temp: " + string_wallace_melting_temp + "\n";
                out += "      SantaLucia Melting Temp: " + string_sl_melting_temp + "\n";
                out += " GC Content of Binding Region: " + string_GC_content + + "\n";
                out += "                     GC Clamp: " + string_GC_clamp + "\n";
                out += "             Homodimerization: " + string_homo_dimerizes + "\n";
                out += "       Homopolymer Run Length: " + string_max_length_homopolymer_run + "\n";
                out += "             Forms Hairpin(s): " + string_hairpin + "\n";
                out += "----------------------------------------------------------------------\n";
                
                //std::cout << out;
                outputFile << out;
            }
            
            /* Print original gene sequence to output .txt file. */
            outputFile << "\n\n"
                       << "|-------------------------|\n"
                       << "| Original Gene Sequence: |\n"
                       << "|-------------------------|\n"
                       << "Length : " << gene_seq.size() << "\n"
                       << gene_seq;
            
            /* Print new mutated gene sequence to output .txt file. */
            outputFile << "\n\n"
                       << "|-------------------------|\n"
                       << "| Mutated Gene Sequence:  |\n"
                       << "|-------------------------|\n"
                       << "Length : " << mutated_gene_seq.size() << "\n"
                       << mutated_gene_seq
                       << "\n";
            
            // Close file
            outputFile.close();

            std::cout << "\nSuccess: A list of viable primers has been output to " << primer_filename.substr(7);
            
        } else {
            std::cerr << "\nError creating primer output file: " << primer_filename.substr(7) << "\n";
            return 1;
        }
        
        
        /* Generate a .fasta file containing the new mutated gene sequence. */
        // Find an available filename.
        std::string mutated_gene_filename = find_available_filename("output/mutated_gene", ".fasta");
        // Create and open a text file for writing in the output folder.
        std::ofstream mutated_gene_output_file(mutated_gene_filename);
        // Check if the file was successfully opened.
        if (mutated_gene_output_file.is_open()) {
            /* Write descriptive fasta header. */
            std::string mutated_gene_header = gene_header;
            mutated_gene_header += " ";
            mutated_gene_header += user_mutation;
            
            // Write mutated sequence to file.
            mutated_gene_output_file << mutated_gene_header << "\n" << mutated_gene_seq;
            
            // Close file.
            mutated_gene_output_file.close();
            
            // Print success message to user.
            std::cout << "\nThe mutated gene sequence has been output to " << mutated_gene_filename.substr(7) << "\n"
                      << "\nMake sure to select your two primers and ensure they do not dimerize by using function (D) of this program.\n\n";
            
            // Success.
            return 0;
        } else {
            std::cerr << "\nError creating mutated gene output file: " << mutated_gene_filename.substr(7) << "\n";
            return 1;
        }
    }
    /*------------------------------- END MUTAGENIC PRIMERS ---------------------------------------- */
    /* --------------------------------------------------------------------------------------------- */

    /* --------------------------------------------------------------------------------------------- */
    
    /* --------------------------------------------------------------------------------------------- */
    /*------------------------------ GENE INSERTION PRIMERS ---------------------------------------- */
    
    else if (user_input == "C" || user_input == "c" || user_input == "insertion") {
        /* C.) GENERATE PRIMERS FOR INSERTION OF A GENE INTO A VECTOR. */
        std::cout << "\n- GENERATING PRIMERS FOR INSERTION OF A GENE INTO A VECTOR -";
        
        // Read in vector and gene FASTA
        std::cout << "\nEnter vector FASTA filename: ";
        std::getline(std::cin, vector_filename);
        
        // Read file and generate vector SequenceRecord object. (Not supported yet)
        // SequenceRecord vector_sr = SequenceRecord::load_from_file(relative_path + vector_filename);
        
        if (!read_fasta(relative_path + vector_filename, vector_seq, vector_header)) {
            std::cerr << "Error reading vector file.\n";
            return 1;
        }
        std::string vector_name = extract_vector_name(vector_header);

        std::cout << "Enter gene FASTA filename: ";
        std::getline(std::cin, gene_filename);
        
        if (!read_fasta(relative_path + gene_filename, gene_seq, gene_header)) {
            std::cerr << "Error reading gene file.\n";
            return 1;
        }
        
        std::vector<RestrictionEnzyme> enzymes = read_RE_json(RE_sites_file_path);
        
        std::cout << "\nLoaded " << enzymes.size() << " restriction enzymes from " << RE_sites_file_path.substr(5) << "\n";
        
        
        /* Organize enzymes that cut the vector and do NOT cut the gene into new C++ vectors. */
        // Restriction enzymes that cut the vector ONCE
        std::vector<RestrictionEnzyme> vector_cutters;
        // Restriction enzymes that do NOT cut the gene
        std::vector<RestrictionEnzyme> gene_non_cutters;
        // Iterate over the list of restriction enzymes
        for (RestrictionEnzyme& e : enzymes) {
            size_t vector_forward_RE_site_index = nucleotide::find_DNA_seq(e.forward_sequence, vector_seq);
            size_t vector_revcomp_RE_site_index = nucleotide::find_DNA_seq(e.reverse_complement, vector_seq);
            size_t gene_forward_RE_site_index = nucleotide::find_DNA_seq(e.forward_sequence, gene_seq);
            size_t gene_revcomp_RE_site_index = nucleotide::find_DNA_seq(e.reverse_complement, gene_seq);
            
            // If forward sequence of e is found in vector ONCE, update e
            // (if forward sequence of e is found MORE than once in vector, e.forward_cuts_vector is left as false)
            if (vector_forward_RE_site_index != std::string::npos && !found_again(vector_forward_RE_site_index, e.forward_sequence, vector_seq)) {
                e.forward_cuts_vector = true;
                e.forward_RE_site_index = vector_forward_RE_site_index;
            }
            
            // If reverse complement of e is found in vector ONCE, update e
            // (if reverse complement of e is found MORE than once in vector, e.revcomp_cuts_vector is left as false)
            if (vector_revcomp_RE_site_index != std::string::npos && !found_again(vector_revcomp_RE_site_index, e.reverse_complement, vector_seq)) {
                e.revcomp_cuts_vector = true;
                e.revcomp_RE_site_index = vector_revcomp_RE_site_index;
            }
            
            // If forward sequence of e is NOT found in gene, update e
            if (gene_forward_RE_site_index == std::string::npos) {
                e.forward_cuts_gene = false;
                gene_non_cutters.push_back(e);
            } 
            
            // If reverse complement of e is NOT found in gene, update e
            if (gene_revcomp_RE_site_index == std::string::npos) {
                e.revcomp_cuts_gene = false;
                gene_non_cutters.push_back(e);
            }
            
            
            // If e cuts vector ONCE, add updated e to vector_cutters
            // 1. Check to see if e cuts vector at all
            if (e.forward_cuts_vector || e.revcomp_cuts_vector) {
                // 2. Check for edge case when e.forward_sequence cuts vector AND e.reverse_complement cuts vector
                if (e.forward_cuts_vector && e.revcomp_cuts_vector) {
                    // 3. If they cut at the same place, they are palindromic
                    if (e.forward_RE_site_index == e.revcomp_RE_site_index) {
                        // Everything is satisfied; add updated e to vector_cutters
                        vector_cutters.push_back(e);
                    }
                }
            }
            
            // If e does NOT cut gene, add updated e to gene_non_cutters
            if (!e.forward_cuts_gene && !e.revcomp_cuts_gene) {
                gene_non_cutters.push_back(e);
            }
        }
        
        // approved_enzymes: restriction enzymes that cut the vector ONCE and do NOT cut the gene
        std::vector<RestrictionEnzyme> approved_enzymes;
        // Generate approved_enzymes
        for (RestrictionEnzyme cutter : vector_cutters) {
            if (!cutter.forward_cuts_gene && !cutter.revcomp_cuts_gene) {
                approved_enzymes.push_back(cutter);
                
            }
        }
        
        
        // Generate .txt file containing the approved restriction enzymes
        // i.e. enzymes that cut the vector ONCE but do NOT cut the gene
        std::string gene_filename_sans_fasta = gene_filename.substr(0, gene_filename.find(".fasta"));
        std::string vector_filename_sans_fasta = vector_filename.substr(0, vector_filename.find(".fasta"));
        std::string approved_RE_filename = find_available_filename("output/viable_restriction_enzymes", ".txt");
        // Create and open a text file for writing
        std::ofstream outputFile(approved_RE_filename);
        // Check if the file was successfully opened
        if (outputFile.is_open()) {
            // Write data to file
            outputFile << " - VIABLE RESTRICTION ENZYMES FOR GENE INSERTION -\n\n"
                       << "(i.e. restriction enzymes that cut the vector ONCE but do NOT cut the gene)\n"
                       << "Gene file: " << gene_filename << "\n"
                       << "Vector file: " << vector_filename << "\n\n";
            
            for (RestrictionEnzyme re : approved_enzymes) {
                std::string out = re.to_string();
                // std::cout << out;
                outputFile << out;
            }
            
            // Close  file
            outputFile.close();

            std::cout << "\nA list of viable restriction enzymes has been output to: " << approved_RE_filename.substr(7);
        } else {
            std::cerr << "\nError opening restriction enzyme output file: " << approved_RE_filename;
            return 1;
        }
        
        // Write linear cut map to output file.
        write_linear_cut_map(vector_name, vector_seq, approved_enzymes);
        
        
        // Create .svg circular map of restriction enzyme cut sites.
        write_svg_vector_map(vector_name, vector_seq, approved_enzymes);
        
        
        /* Query user which 2 restriction enzymes they would like to use. */
        std::string user_enzyme_name1;
        RestrictionEnzyme user_enzyme1;
        std::string user_enzyme_name2;
        RestrictionEnzyme user_enzyme2;
        bool valid_enzyme = false;
        while (!valid_enzyme) {
            std::cout << "\n\nEnter the name of the first (5') (upstream) restriction enzyme you would like to use: ";
            std::getline(std::cin, user_enzyme_name1);
            for (RestrictionEnzyme e : approved_enzymes) {
                if (e.name == user_enzyme_name1) {
                    valid_enzyme = true;
                    user_enzyme1 = e;  // shallow copy is fine.
                }
            }
            if (!valid_enzyme) {
                std::cout << "That enzyme was not recognized. Make sure to copy paste the entire name including any isoschizomers after the space. Please try again.\n";
            }
        }
        bool valid_enzyme2 = false;
        while (!valid_enzyme2) {
            std::cout << "\nEnter the name of the second (3') (downstream) restriction enzyme you would like to use: ";
            std::getline(std::cin, user_enzyme_name2);
            for (RestrictionEnzyme e : approved_enzymes) {
                if (e.name == user_enzyme_name2) {
                    valid_enzyme2 = true;
                    user_enzyme2 = e;  // shallow copy is fine.
                }
            }
            if (!valid_enzyme2) {
                std::cout << "That enzyme was not recognized. Make sure to copy paste the entire name including any isoschizomers after the space. Please try again.\n";
            }
        }
        
        /* Query user primer specifications.
           A.) Forward primer specifications:
            1. Add N-terminus tag?
                [YES]: Select tag: His6, FLAG, HA, Myc, other.
            2. Add N-terminus protease?
                [YES]: Select protease: TEV, PreScission, FXa, Enterokinase, other.
            3. Add start codon?
                [NO]: a. Enter the position of the first n.t. of the vector start codon.
                      b. Add flexible linker?
                          [YES]: Enter linker sequence.
                      c. (IFF R.F. OUT OF PHASE) Reading frame is 1/2 n.t. out of phase. Spacers ok?
                          [NO]: Please choose a new 5' restriction enzyme.
         
           B.) Reverse primer specifications:
            1. Add C-terminus tag?
                [YES]: Select tag: His6, FLAG, HA, Myc, other. //Remember revcomp
            2. Add C-terminus protease?
                [YES]: Select protease: TEV, PreScission, FXa, Enterokinase, other.
            3. Add stop codon?
                [NO]: a. Enter the position of the first n.t. of the subsequent vector gene.
                      b. Add flexible linker?
                          [YES]: Enter linker sequence.
                      c. (IFF R.F. OUT OF PHASE) Reading frame is 1/2 n.t. out of phase. Spacers ok?
                          [NO]: Please choose a new 3' restriction enzyme.
         */
        struct PrimerSpecs {
            // ------------------------------------ Forward
            bool want_n_tag = false;
            std::string n_tag_seq = "";
            
            bool want_n_protease = false;
            std::string n_protease_seq = "";
            
            bool want_start = false;
            std::string start_seq = "ATG";
            std::size_t index_vector_atg_0_based = std::string::npos;
            
            bool want_n_linker = false;
            std::string n_linker_seq = "";
            
            bool want_n_spacer;
            std::string n_spacer_seq = "";
            // ------------------------------------
            // ------------------------------------ Reverse
            bool want_c_tag = false;
            std::string c_tag_seq = "";
            
            bool want_c_protease = false;
            std::string c_protease_seq = "";
            
            bool want_stop = false;
            std::string stop_seq = "TGA";
            std::size_t index_vector_next_gene_0_based = std::string::npos;
            
            bool want_c_linker = false;
            std::string c_linker_seq = "";
            
            bool want_c_spacer = false;
            std::string c_spacer_seq = "";
            // ------------------------------------
        };
        
        PrimerSpecs primer_specs;
        
        // Forward primer specs:
        /* -------------------------- N-terminus tag -------------------------- */
        output = "Add N-terminus tag? (y/n)";
        valid_inputs = {"y", "n", "q"};
        error_message = "Unrecognized option.";
        user_input = query_user(output, valid_inputs, error_message);
        if (user_input == "q") { // escape sequence.
            std::cout << "\nQuitting program.\n";
            return 1;
            
        } else if (user_input == "y") { // yes n-tag.
            primer_specs.want_n_tag = true;
            
            output = "Select tag. (His6/FLAG/HA/Myc/Gst/other)";
            valid_inputs = {"His6", "FLAG", "HA", "Myc", "other", "x"};
            error_message = "Unrecognized option.";
            user_input = query_user(output, valid_inputs, error_message);
            if (user_input == "q") { // escape sequence
                std::cout << "\nQuitting program.\n";
                return 1;
                
            } else if (user_input == "other") { // custom tag
                output = "Enter tag DNA sequence.";
                valid_inputs = {};
                error_message = "";
                user_input = query_user(output, valid_inputs, error_message);
            }
            primer_specs.n_tag_seq = generate_tag_seq(user_input);
            
        } else if (user_input == "n") { // no n-tag.
            primer_specs.want_n_tag = false;

        }
        
        /* -------------------------- N-terminus protease -------------------------- */
        output = "Add N-terminus protease? (y/n)";
        valid_inputs = {"y", "n", "q"};
        error_message = "Unrecognized option.";
        user_input = query_user(output, valid_inputs, error_message);
        if (user_input == "y") { // yes n-protease
            primer_specs.want_n_protease = true;
            
            /* -------------------------- Protease selection -------------------------- */
            output = "Select protease. (TEV, PreScission, FXa, Enterokinase, other)";
            valid_inputs = {"TEV", "PreScission", "FXa", "Enterokinase", "other", "q"};
            error_message = "Unrecognized option.";
            user_input = query_user(output, valid_inputs, error_message);
            if (user_input == "q") { // escape sequence.
                std::cout << "\nQuitting program.\n";
                return 1;
            }
            if (user_input == "other") { // custom protease sequence
                output = "Enter protease recognition sequence";
                valid_inputs = {};
                error_message = "";
                user_input = query_user(output, valid_inputs, error_message);
            }
            primer_specs.n_protease_seq = generate_protease_seq(user_input);
            /* --------------------------------------------------------------------- */
            
        } else if (user_input == "n") { // no n-protease.
            primer_specs.want_n_protease = false;

        } else if (user_input == "q") { // escape sequence.
            std::cout << "\nQuitting program.\n";
            return 1;
        }

        /* -------------------------- Start codon -------------------------- */
        output = "Add start codon? (y/n)";
        valid_inputs = {"y", "n", "q"};
        error_message = "Unrecognized option.";
        user_input = query_user(output, valid_inputs, error_message);
        if (user_input == "q") { // escape sequence
            std::cout << "\nQuitting program.\n";
            return 1;
        } else if (user_input == "y") { // yes start codon.
            primer_specs.want_start = true;
            
        } else if (user_input == "n") { // no start codon.
            primer_specs.want_start = false;
            primer_specs.start_seq = "";
            
            /* ------------------- Position of vector start codon ------------------ */
            output = "Enter the position of the first n.t. of the vector start codon";
            valid_inputs = {};
            error_message = "";
            user_input = query_user(output, valid_inputs, error_message);
            if (user_input == "q") { // escape sequence.
                std::cout << "\nQuitting program.\n";
                return 1;
            }
            // Convert position of start from string to
            //  size_t and convert to 0-base index. */
            std::stringstream ss(user_input);
            size_t user_input_size_t;
            ss >> user_input_size_t;
            primer_specs.index_vector_atg_0_based = user_input_size_t - 1;
            
            /* -------------------------- Flexible n-linker ------------------------ */
            output = "Add N-terminus flexible linker? (y/n)";
            valid_inputs = {"y", "n", "q"};
            error_message = "Unrecognized option.";
            user_input = query_user(output, valid_inputs, error_message);
            if (user_input == "q") { // escape sequence.
                std::cout << "\nQuitting program.\n";
                return 1;
            }
            if (user_input == "y") { // yes flexible linker.
                primer_specs.want_n_linker = true;
                /* ------------------- Choose sequence of flexible n-linker  ------------------- */
                output = "Enter N-terminus linker sequence";
                valid_inputs = {};
                error_message = "";
                user_input = query_user(output, valid_inputs, error_message);
                // Check for escape sequence.
                if (user_input == "q") {
                    std::cout << "\nQuitting program.\n";
                    return 1;
                }
                // Assign linker sequence.
                primer_specs.n_linker_seq = user_input;
                /* --------------------------------------------------------------------------- */
 
            } else if (user_input == "n") { // no flexible n-linker.
                primer_specs.want_n_linker = false;
            }
            
            /* ------------------------------- Spacers ----------------------------- */
            // Compute reading frame phase of first gene nucleotide.
            // Use the index of the vector start codon as r.f. anchor.

            //   Vector start codon -...- RE sequence - Flexible linker - [Spacer?] - Gene
            //   |---------------------|
            size_t distance_start_to_RE_seq = user_enzyme1.forward_RE_site_index - primer_specs.index_vector_atg_0_based;
            
            //   Vector start codon -...- RE sequence - Flexible linker - [Spacer?] - Gene
            //                           |-----------|
            size_t RE_seq_length = user_enzyme1.forward_sequence.size();
            
            //   Vector start codon -...- RE sequence - Flexible linker - [Spacer?] - Gene
            //                                         |---------------|
            size_t linker_length = primer_specs.n_linker_seq.size();
            
            //   Vector start codon -...- RE sequence - Flexible linker - [Spacer?] - Gene
            //   |-----------------------------------------------------|
            size_t distance_start_to_gene = distance_start_to_RE_seq + RE_seq_length + linker_length;
            
            size_t phase_gene = distance_start_to_gene % 3;
            if (phase_gene > 0) {
                size_t length_spacer_needed = 3 - phase_gene;

                output = "Gene is out of phase. Add " + std::to_string(length_spacer_needed) + " nucleotide N-terminus spacer? (y/n)";
                valid_inputs = {"y", "n", "q"};
                error_message = "Unrecognized option.";
                user_input = query_user(output, valid_inputs, error_message);
                if (user_input == "q") { // escape sequence.
                    std::cout << "\nQuitting program.\n";
                    return 1;
                } else if (user_input == "n") { // no add n-spacer
                    primer_specs.want_n_spacer = false;
                    std::cout << "\nPlease choose a new 5' restriction enzyme.\n";
                    return 1;
                } else if (user_input == "y") { // yes add n-spacer
                    primer_specs.want_n_spacer = true;
                    /* ------------------------------- Determine spacers ----------------------------- */
                    output = "Enter N-terminus spacer sequence (" + std::to_string(length_spacer_needed) + " nucleotides)";
                    valid_inputs = {"#length:"+std::to_string(length_spacer_needed), "q"};
                    error_message = "Wrong length.";
                    user_input = query_user(output, valid_inputs, error_message);
                    if (user_input == "q") { // escape sequence.
                        std::cout << "\nQuitting program.\n";
                        return 1;
                    } else { // user input is spacer sequence.
                        primer_specs.n_spacer_seq = nucleotide::normalize(user_input);
                    }
                    
                }
            }
        } // End of "don't add start codon" logic.
        // End of forward primer specs.
        
        
        // Reverse primer specs:
        /* -------------------------- C-terminus tag -------------------------- */
        output = "Add C-terminus tag? (y/n)";
        valid_inputs = {"y", "n", "q"};
        error_message = "Unrecognized option.";
        user_input = query_user(output, valid_inputs, error_message);
        if (user_input == "q") { // escape sequence.
            std::cout << "\nQuitting program.\n";
            return 1;
            
        } else if (user_input == "y") { // yes c-tag.
            primer_specs.want_c_tag = true;
            
            output = "Select tag. (His6/FLAG/HA/Myc/Gst/other)";
            valid_inputs = {"His6", "FLAG", "HA", "Myc", "other", "q"};
            error_message = "Unrecognized option.";
            user_input = query_user(output, valid_inputs, error_message);
            if (user_input == "q") { // escape sequence
                std::cout << "\nQuitting program.\n";
                return 1;
                
            } else if (user_input == "other") { // custom tag
                output = "Enter tag DNA sequence.";
                valid_inputs = {};
                error_message = "";
                user_input = query_user(output, valid_inputs, error_message);
            }
            primer_specs.c_tag_seq = generate_tag_seq(user_input);
            
        } else if (user_input == "n") { // no c-tag.
            primer_specs.want_c_tag = false;

        }
        
        /* -------------------------- C-terminus protease -------------------------- */
        output = "Add C-terminus protease? (y/n)";
        valid_inputs = {"y", "n", "x"};
        error_message = "Unrecognized option.";
        user_input = query_user(output, valid_inputs, error_message);
        if (user_input == "y") { // yes c-protease
            primer_specs.want_c_protease = true;
            
            /* -------------------------- Protease selection -------------------------- */
            output = "Select protease. (TEV, PreScission, FXa, Enterokinase, other)";
            valid_inputs = {"TEV", "PreScission", "FXa", "Enterokinase", "other", "q"};
            error_message = "Unrecognized option.";
            user_input = query_user(output, valid_inputs, error_message);
            if (user_input == "q") { // escape sequence.
                std::cout << "\nQuitting program.\n";
                return 1;
            }
            if (user_input == "other") { // custom protease sequence
                output = "Enter protease recognition sequence";
                valid_inputs = {};
                error_message = "";
                user_input = query_user(output, valid_inputs, error_message);
            }
            primer_specs.c_protease_seq = generate_protease_seq(user_input);
            /* --------------------------------------------------------------------- */
            
        } else if (user_input == "n") { // no c-protease.
            primer_specs.want_c_protease = false;

        } else if (user_input == "q") { // escape sequence.
            std::cout << "\nQuitting program.\n";
            return 1;
        }

        /* -------------------------- Stop codon -------------------------- */
        output = "Add stop codon? (y/n)";
        valid_inputs = {"y", "n", "q"};
        error_message = "Unrecognized option.";
        user_input = query_user(output, valid_inputs, error_message);
        if (user_input == "q") { // escape sequence
            std::cout << "\nQuitting program.\n";
            return 1;
        } else if (user_input == "y") { // yes stop codon.
            primer_specs.want_stop = true;
            
        } else if (user_input == "n") { // no stop codon.
            primer_specs.want_stop = false;
            primer_specs.stop_seq = "";
            
            /* ------------------- Position of next gene in vector ------------------ */
            output = "Enter the position of the beginning of the next gene in the vector";
            valid_inputs = {};
            error_message = "";
            user_input = query_user(output, valid_inputs, error_message);
            if (user_input == "q") { // escape sequence.
                std::cout << "\nQuitting program.\n";
                return 1;
            }
            // Convert position of next gene from string to
            //  size_t and convert to 0-base index. */
            std::stringstream ss(user_input);
            size_t user_input_size_t;
            ss >> user_input_size_t;
            primer_specs.index_vector_next_gene_0_based = user_input_size_t - 1;
            
            /* -------------------------- Flexible c-linker ------------------------ */
            output = "Add C-terminus flexible linker? (y/n)";
            valid_inputs = {"y", "n", "q"};
            error_message = "Unrecognized option.";
            user_input = query_user(output, valid_inputs, error_message);
            if (user_input == "q") { // escape sequence.
                std::cout << "\nQuitting program.\n";
                return 1;
            }
            if (user_input == "y") { // yes c-flexible linker.
                primer_specs.want_c_linker = true;
                /* ------------------- Choose sequence of flexible c-linker  ------------------- */
                output = "Enter C-terminus linker sequence";
                valid_inputs = {};
                error_message = "";
                user_input = query_user(output, valid_inputs, error_message);
                // Check for escape sequence.
                if (user_input == "q") {
                    std::cout << "\nQuitting program.\n";
                    return 1;
                }
                // Assign linker sequence.
                primer_specs.c_linker_seq = user_input;
                /* --------------------------------------------------------------------------- */
 
            } else if (user_input == "n") { // no flexible c-linker.
                primer_specs.want_c_linker = false;
            }
            
            /* ------------------------------- Spacers ----------------------------- */
            // Compute reading frame phase of gene to next gene.
            // Assume gene is placed perfectly in frame.

            //   Gene - Protease - Tag - Linker - [Spacers?] - RE2 - ... - Next gene
            //   |---|
            //   gene_seq.size()
            //         |--------|
            //         primer_specs.c_protease_seq.size()
            //                    |----|
            //                    primer_specs.c_tag_seq.size()
            //                          |------|
            //                          primer_specs.c_linker_seq.size()
            //                                                |---|
            //                                                user_enzyme2.forward_sequence.size()

            //   Gene - Protease - Tag - Linker - [Spacers?] - RE2 - ... - Next gene
            //   |------------------------------------------------|
            size_t gene_to_RE2_length = gene_seq.size()
                                      + primer_specs.c_protease_seq.size()
                                      + primer_specs.c_tag_seq.size()
                                      + primer_specs.c_linker_seq.size()
                                      + user_enzyme2.forward_sequence.size();
            
            //   Gene - Protease - Tag - Linker - [Spacers?] - RE2 - ... - Next gene
            //                                                      |---|
            size_t distance_RE2_next_gene = primer_specs.index_vector_next_gene_0_based
                                          - user_enzyme2.forward_RE_site_index
                                          - user_enzyme2.forward_sequence.size();
            
            //   Gene - Protease - Tag - Linker - [Spacers?] - RE2 - ... - Next gene
            //   |------------------------------------------------------|
            size_t gene_to_next_gene_length = gene_to_RE2_length + distance_RE2_next_gene;

            size_t phase_next_gene = gene_to_next_gene_length % 3;
            if (phase_next_gene > 0) {
                size_t length_c_spacer_needed = 3 - phase_next_gene;

                output = "Next gene is out of phase. Add " + std::to_string(length_c_spacer_needed) + " nucleotide C-terminus spacer? (y/n)";
                valid_inputs = {"y", "n", "q"};
                error_message = "Unrecognized option.";
                user_input = query_user(output, valid_inputs, error_message);
                if (user_input == "q") { // escape sequence.
                    std::cout << "\nQuitting program.\n";
                    return 1;
                } else if (user_input == "n") { // no add c-spacer
                    primer_specs.want_c_spacer = false;
                    std::cout << "\nPlease choose a new 3' restriction enzyme.\n";
                    return 1;
                } else if (user_input == "y") { // yes add c-spacer
                    primer_specs.want_c_spacer = true;
                    /* ------------------------------- Determine spacers ----------------------------- */
                    output = "Enter C-terminus spacer sequence (" + std::to_string(length_c_spacer_needed) + " nucleotides)";
                    valid_inputs = {"#length:"+std::to_string(length_c_spacer_needed), "q"};
                    error_message = "Wrong length.";
                    user_input = query_user(output, valid_inputs, error_message);
                    if (user_input == "q") { // escape sequence.
                        std::cout << "\nQuitting program.\n";
                        return 1;
                    } else { // user input is spacer sequence.
                        primer_specs.c_spacer_seq = nucleotide::normalize(user_input);
                    }
                }
            }
        } // End of "don't add stop codon" logic.
        // End of reverse primer specs
        
        
        /* Create insertion primers based on selected restriction enzymes and user input.
           Primer architecture:
         Forward:
         5' - [clamp] [RE1] [spacers?] [N-tag?] [protease?] [linker/start?] [annealing] - 3'
         Reverse: NOTE "revcomp" is abbreviated "RC"
         5' - RC[clamp] RC[RE2] RC[spacers?] RC[linker/stop?] RC[C-tag?] RC[protease?] [annealing] - 3'
         */
        // Default clamp to provide space for restriction enzyme to successfully cut.
        const std::string CLAMP_SEQ = "GCGTAC";
        
        // Non annealing part of the forward primer.
        // [clamp] [RE1] [spacers?] [N-tag?] [protease?] [linker/start?]
        std::string f_non_annealing_seq = CLAMP_SEQ
                                        + user_enzyme1.forward_sequence
                                        + primer_specs.n_spacer_seq
                                        + primer_specs.n_tag_seq
                                        + primer_specs.n_protease_seq
                                        + primer_specs.start_seq
                                        + primer_specs.n_linker_seq;
        
        // Non annealing part of the reverse primer.
        // RC[clamp] RC[RE2] RC[spacers?] RC[linker/stop?] RC[C-tag?] RC[protease?]
        std::string r_non_annealing_seq = nucleotide::revcomp_of(CLAMP_SEQ)
                                        + nucleotide::revcomp_of(user_enzyme2.forward_sequence)
                                        + nucleotide::revcomp_of(primer_specs.c_spacer_seq)
                                        + nucleotide::revcomp_of(primer_specs.c_linker_seq)
                                        + nucleotide::revcomp_of(primer_specs.stop_seq)
                                        + nucleotide::revcomp_of(primer_specs.c_tag_seq)
                                        + nucleotide::revcomp_of(primer_specs.c_protease_seq);
        
        // Make beginning of forward primer notation string.
        std::string beg_f_notation = "|Clamp|";
        beg_f_notation += "|RE1" + repeat(" ", user_enzyme1.forward_sequence.size() - 5) + "|";
        if (primer_specs.want_n_spacer)
            beg_f_notation += repeat("_", primer_specs.n_spacer_seq.size());
        if (primer_specs.want_n_tag)
            beg_f_notation += "|Tag" + repeat(" ", primer_specs.n_tag_seq.size() - 5) + "|";
        if (primer_specs.want_n_protease)
            beg_f_notation += "|Protease" + repeat(" ", primer_specs.n_protease_seq.size() - 10) + "|";
        if (primer_specs.want_start)
            beg_f_notation += "|>|";
        if (primer_specs.want_n_linker) {
            if (primer_specs.n_linker_seq.size() >= 8)
                beg_f_notation += "|Linker" + repeat(" ", primer_specs.n_linker_seq.size() - 8) + "|";
            if (primer_specs.n_linker_seq.size() >= 3 && primer_specs.n_linker_seq.size() < 8)
                beg_f_notation += "|L" + repeat(" ", primer_specs.n_linker_seq.size() - 3) + "|";
            if (primer_specs.n_linker_seq.size() < 3)
                beg_f_notation += repeat("_", primer_specs.n_linker_seq.size());
        }
        beg_f_notation += "|Annealing";
        
        // Make beginning of reverse primer notation string.
        // The string "|Annealing ...... |" will be prepended later when we know the length of the annealing region.
        std::string beg_r_notation = "|Clamp|";
        beg_r_notation += "|RE2" + repeat(" ", user_enzyme2.forward_sequence.size() - 5) + "|";
        if (primer_specs.want_c_spacer)
            beg_r_notation += repeat("_", primer_specs.c_spacer_seq.size());
        if (primer_specs.want_c_linker) {
            if (primer_specs.c_linker_seq.size() >= 8)
                beg_r_notation += "|Linker" + repeat(" ", primer_specs.c_linker_seq.size() - 8) + "|";
            if (primer_specs.c_linker_seq.size() >= 3 && primer_specs.c_linker_seq.size() < 8)
                beg_r_notation += "|L" + repeat(" ", primer_specs.c_linker_seq.size() - 3) + "|";
            if (primer_specs.c_linker_seq.size() < 3)
                beg_r_notation += repeat("_", primer_specs.c_linker_seq.size());
        }
        if (primer_specs.want_stop)
            beg_r_notation += "|*|";
        if (primer_specs.want_c_tag)
            beg_r_notation += "|Tag" + repeat(" ", primer_specs.c_tag_seq.size() - 5) + "|";
        if (primer_specs.want_c_protease)
            beg_r_notation += "|Protease" + repeat(" ", primer_specs.c_protease_seq.size() - 10) + "|";
        beg_r_notation += "|Annealing";


        // Vectors of Primer objects representing annealing parts of each primer.
        std::vector<Primer> forward_annealing_portions; // Will be ordered based on length
        std::vector<Primer> reverse_annealing_portions; // "
                
        // Generate list of viable FORWARD annealing primer portions (forward_annealing_portions)
        for (size_t i = 18; i <= 25; i++) {
            // Forward primers bind to the reverse strand
            std::string primer_seq = gene_seq.substr(0, i);
            // Create a potential primer
            Primer x = Primer(primer_seq);
            // Determine if x satisfies all internal rules
            if (x.is_good()) {
                // Ensure x doesn't bind off-target
                if (!x.binds_off_target(gene_seq)) {
                    // All internal rules met + Primer does not bind off-target; x added to list of viable FORWARD annealing portions
                    forward_annealing_portions.push_back(x);
                }
            }
        }
        // Generate list of viable REVERSE annealing primer portions (reverse_annealing_portions)
        for (size_t i = 18; i <= 25; i++) {
            // Note: Reverse primers bind to the forward strand.
            std::string forward_end = gene_seq.substr(gene_seq.size() - i);
            // Reverse complement of the end of the forward seqeunce of the gene is the reverse primer.
            std::string reverse_primer_seq = nucleotide::revcomp_of(forward_end);
            // Create a potential primer
            Primer x = Primer(reverse_primer_seq);
            // Determine if x satisfies all internal rules
            if (x.is_good()) {
                // Ensure x doesn't bind off-target
                if (!x.binds_off_target(gene_seq)) {
                    // All internal rules met + Primer does not bind off-target; x added to list of viable REVERSE annealing portions
                    reverse_annealing_portions.push_back(x);
                }
            }
        }

        /* Combine non-annealing portion with annealing portions and add to lists of full FORWARD primer 
            sequences.
         */
        std::vector<std::string> forward_full_primer_seqs;
        for (Primer annealing_portion: forward_annealing_portions) {
            std::string full_primer_seq = f_non_annealing_seq
                                          + annealing_portion.seq;
            forward_full_primer_seqs.push_back(full_primer_seq);
        }
        /* Combine non-annealing portion with annealing portions, generate revcomp, and add revcomp to
            lists of full REVERSE primer sequences.
         */
        std::vector<std::string> reverse_full_primer_seqs;
        for (Primer annealing_portion: reverse_annealing_portions) {
            // RC[clamp] RC[RE2] RC[spacers?] RC[linker/stop?] RC[C-tag?] RC[protease?]
            std::string full_primer_seq = r_non_annealing_seq
                                          + annealing_portion.seq;
            reverse_full_primer_seqs.push_back(full_primer_seq);
        }
        
        
        
        // Generate .txt file containing the approved forward and reverse primers for gene insertion.
        // Find an available filename.
        std::string gene_ins_primer_filename = find_available_filename("output/gene_insertion_primers", ".txt");
        // Create and open a text file for writing in the output folder.
        std::ofstream outputFile2(gene_ins_primer_filename);
        // Check if the file was successfully opened.
        if (outputFile2.is_open()) {
            // Write data to file.
            outputFile2 << " - VIABLE PRIMERS FOR GENE INSERTION INTO VECTOR -\n\n"
                       << "Gene source file: " << gene_filename << "\n"
                       << "Vector source file: " << vector_filename << "\n\n"
                       << "5' Restriction enzyme: " << user_enzyme1.name
                       << " (" << user_enzyme1.notated_sequence << ")" << "\n"
                       << "3' Restriction enzyme Name: " << user_enzyme2.name
                       << " (" << user_enzyme2.notated_sequence << ")" << "\n";


            // FORWARD primers write to output file.
            std::string forward_primer_title = "\n"
                                     "|--------------------------------------------|\n"
                                     "| Viable FORWARD Primers for gene insertion: |\n"
                                     "|--------------------------------------------|\n";
            outputFile2 << forward_primer_title;
            for (size_t i = 0; i < forward_full_primer_seqs.size(); i++) {
                
                Primer* p = &forward_annealing_portions[i];
                
                std::string string_seq = "5'- " + forward_full_primer_seqs[i] + " -3'";
                std::string string_length = std::to_string(forward_full_primer_seqs[i].size()) + " n.t.";
                std::string string_wallace_melting_temp = std::to_string(p->wallace_melting_temp) + "°C";
                std::string string_sl_melting_temp = std::to_string(p->santalucia_melting_temp).substr(0, 4) + "°C";
                std::string string_GC_content = std::to_string(p->GC_content*100.).substr(0, 4) + "%";
                std::string string_GC_clamp = p->GC_clamp ? "YES" : "NO";
                std::string string_homo_dimerizes = p->homo_dimerizes ? "YES" : "NO";
                std::string string_max_length_homopolymer_run = std::to_string(p->max_length_homopolymer_run) + " n.t.";
                std::string string_hairpin = p->hairpin ? "YES" : "NO";
               
                // Append spaces and "|" to f_notation.
                std::string f_notation = beg_f_notation + repeat(" ", p->seq.size() - 10) + "|";
                
                std::string out = "";
                out += "----------------------------------------------------------------------\n";
                out += "              Primer Sequence: " + string_seq + "\n";
                out += "                                  " + f_notation + "\n";
                out += "                       Length: " + string_length + "\n";
                out += "         Wallace Melting Temp: " + string_wallace_melting_temp + "\n";
                out += "      SantaLucia Melting Temp: " + string_sl_melting_temp + "\n";
                out += " GC Content of Binding Region: " + string_GC_content + + "\n";
                out += "                     GC Clamp: " + string_GC_clamp + "\n";
                out += "             Homodimerization: " + string_homo_dimerizes + "\n";
                out += "       Homopolymer Run Length: " + string_max_length_homopolymer_run + "\n";
                out += "             Forms Hairpin(s): " + string_hairpin + "\n";
                out += "----------------------------------------------------------------------\n";
                
                //std::cout << out;
                outputFile2 << out;
            }
            
            // REVERSE primers write to output file
            std::string reverse_primer_title = "\n\n"
                                     "|--------------------------------------------|\n"
                                     "| Viable REVERSE Primers for gene insertion: |\n"
                                     "|--------------------------------------------|\n";
            outputFile2 << reverse_primer_title;
            for (size_t i = 0; i < reverse_full_primer_seqs.size(); i++) {
                
                Primer* r = &reverse_annealing_portions[i];
                
                std::string string_seq = "5'- " + reverse_full_primer_seqs[i] + " -3'";
                std::string string_length = std::to_string(reverse_full_primer_seqs[i].size()) + " n.t.";
                std::string string_wallace_melting_temp = std::to_string(r->wallace_melting_temp) + "°C";
                std::string string_sl_melting_temp = std::to_string(r->santalucia_melting_temp).substr(0, 4) + "°C";
                std::string string_GC_content = std::to_string(r->GC_content*100.).substr(0, 4) + "%";
                std::string string_GC_clamp = r->GC_clamp ? "YES" : "NO";
                std::string string_homo_dimerizes = r->homo_dimerizes ? "YES" : "NO";
                std::string string_max_length_homopolymer_run = std::to_string(r->max_length_homopolymer_run) + " n.t.";
                std::string string_hairpin = r->hairpin ? "YES" : "NO";
                
                // Append spaces and "|" to f_notation.
                std::string r_notation = beg_r_notation + repeat(" ", r->seq.size() - 10) + "|";

                std::string out = "";
                out += "----------------------------------------------------------------------\n";
                out += "              Primer Sequence: " + string_seq + "\n";
                out += "                                  " + r_notation + "\n";
                out += "                       Length: " + string_length + "\n";
                out += "         Wallace Melting Temp: " + string_wallace_melting_temp + "\n";
                out += "      SantaLucia Melting Temp: " + string_sl_melting_temp + "\n";
                out += " GC Content of Binding Region: " + string_GC_content + + "\n";
                out += "                     GC Clamp: " + string_GC_clamp + "\n";
                out += "             Homodimerization: " + string_homo_dimerizes + "\n";
                out += "       Homopolymer Run Length: " + string_max_length_homopolymer_run + "\n";
                out += "             Forms Hairpin(s): " + string_hairpin + "\n";
                out += "----------------------------------------------------------------------\n";
                
                //std::cout << out;
                outputFile2 << out;
            }
            
            // Close file.
            outputFile2.close();

            // Success.
            std::cout << "\nSuccess: A list of viable gene insertion primers has been output to "
                      << gene_ins_primer_filename.substr(7)
                      << "\nSelect your two primers and ensure they do not dimerize by using function (D) of this program.\n\n";
            return 0;
            
        } else {
            std::cerr << "\nError creating primer output file: " << gene_ins_primer_filename.substr(7) << "\n";
            return 1;
        }
    }
    /*------------------------------ END GENE INSERTION PRIMERS ------------------------------------ */
    /* --------------------------------------------------------------------------------------------- */
    
    /* --------------------------------------------------------------------------------------------- */
    
    /* --------------------------------------------------------------------------------------------- */
    /*--------------------------------- CHECK DIMERIZATION ----------------------------------------- */
    
    else if (user_input == "D" || user_input == "d" || user_input == "dimerization") {
        /* D.) CHECK FOR DIMERIZATION BETWEEN TWO PRIMERS. */
        std::cout << "\n- CHECKING FOR DIMERIZATION BETWEEN TWO PRIMERS -";
        
        // Flag for ending the program only when the user provides two primers that do not dimerize
        bool dimerizes = true;
        while (dimerizes) {
            // Query user which forward primer they would like to use.
            std::string user_forward_primer_seq;
            std::cout << "\nEnter the sequence of the forward primer (5'->3') you would like to use: ";
            std::getline(std::cin, user_forward_primer_seq);
            
            // Check for escape string ("q")
            if (user_forward_primer_seq == "q") {
                std::cout << "\nQuitting program.\n";
                return 1;
            }
            
            // Normalize forward sequence
            user_forward_primer_seq = nucleotide::normalize(user_forward_primer_seq);
            
            // Query user which reverse primer they would like to use.
            std::string user_reverse_primer_seq;
            std::cout << "\nEnter the sequence of the reverse primer (5'->3') you would like to use: ";
            std::getline(std::cin, user_reverse_primer_seq);
            
            // Check for escape string ("q")
            if (user_reverse_primer_seq == "q") {
                std::cout << "\nQuitting program.";
                return 1;
            }
            
            // Normalize forward sequence
            user_reverse_primer_seq = nucleotide::normalize(user_reverse_primer_seq);
            
            // Create forward primer object
            Primer user_forward_primer = Primer(user_forward_primer_seq);

            // Check for dimerization
            // NOTE: dimerizes_with takes string as param
            dimerizes = user_forward_primer.dimerizes_with(user_reverse_primer_seq);
            
            // User output
            if (dimerizes) {
                std::cout << "\nThese two primers dimerize. Try again with another two primers."
                          << "Enter 'q' to exit."
                          << "\n\n----------------------------------------------------------------------\n";
            } else {
                std::cout << "\nSUCCESS. These two primers DO NOT dimerize.\n\n";
            }
        }
        
        return 0;
    }
    /*------------------------------- END CHECK DIMERIZATION --------------------------------------- */
    /* --------------------------------------------------------------------------------------------- */

    /* --------------------------------------------------------------------------------------------- */
    
    /* --------------------------------------------------------------------------------------------- */
    /*---------------------------------- GENERATE REVCOMP ------------------------------------------ */
    
    else if (user_input == "E" || user_input == "e" || user_input == "revcomp") {
        /* E.) GENERATE REVERSE COMPLEMENT OF A SEQUENCE. */
        std::cout << "\n- GENERATING THE REVERSE COMPLEMENT OF A SEQUENCE -";
        
        std::string user_sequence;
        std::cout << "\nEnter sequence (5'->3') to generate the reverse complement of: ";
        std::getline(std::cin, user_sequence);
        
        user_sequence = nucleotide::normalize(user_sequence);
        
        std::string revcomp_user_sequence = nucleotide::revcomp_of(user_sequence);
        std::cout << "\nThe reverse complement of your sequence is: 5'- "
                  << revcomp_user_sequence
                  << " -3'\n\n";
        
        return 0;
    }
    /*-------------------------------- END GENERATE REVCOMP ---------------------------------------- */
    /* --------------------------------------------------------------------------------------------- */
    
    /* --------------------------------------------------------------------------------------------- */

    /* --------------------------------------------------------------------------------------------- */
    /*------------------------------- UNRECOGNIZED OPTION ------------------------------------------ */
    else {
        /* COULD NOT RECOGNIZE AN OPTION */
        std::cout << "Sorry, that option was not recognized. Ending program.\n\n";
        return 1;
    }
    /*----------------------------- END UNRECOGNIZED OPTION ---------------------------------------- */
    /* --------------------------------------------------------------------------------------------- */


    

    
    return 0;
}

/* ------------------------------------------ END MAIN ----------------------------------------- */
/* --------------------------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------------------------- */
/* ------------------------------------- FUNCTION DEFINITIONS ---------------------------------- */

/*
 * Reads a FASTA file and returns a cleaned sequence string.
 * Ignores lines starting with '>'
 * Removes whitespace, numbers, and other non-letters
 * Converts all letters to uppercase
 */
bool read_fasta(const std::string& filename, std::string& sequence, std::string& header) {
    sequence.clear();
    std::ifstream infile(filename);
    if (!infile) return false;

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') { // FASTA header
            header = line;
            continue;
        }

        // Remove spaces, numbers, punctuation; convert to uppercase
        for (char c : line) {
            if (std::isalpha(static_cast<unsigned char>(c))) {
                sequence += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            }
        }
    }

    return true;
}


/*
 * Reads a JSON file of restriction enzyme names (name), recognition sequences (sequence), 
 * and cut sites (cut:top_strand and cut:bottom_strand)
 * Creates a vector of RestrictionEnzyme objects and assigns all
 */
std::vector<RestrictionEnzyme> read_RE_json(const std::string& filename) {
    std::vector<RestrictionEnzyme> enzymes;
    
    // Create an ifstream object to open the file
    std::ifstream ifs(filename);

    // Check if the file was opened successfully
    if (!ifs.is_open()) {
        std::cerr << "Error: Failed to open JSON file '"
                  << filename
                  << "'. System message: "
                  << std::strerror(errno)
                  << std::endl;
        // Return empty list
        return enzymes;
    }
    // Parse the JSON data from the file stream
    nlohmann::json j;
    try {
        ifs >> j; // Read and parse the JSON directly into the nlohmann::json object
    } catch (const nlohmann::json::parse_error& e) {
        std::cerr << "JSON parse error: " << e.what() << std::endl;
        // Return empty list
        return enzymes;
    }

    for (const auto& enz : j) {
        RestrictionEnzyme e;
        e.name = enz.value("name", "");
        e.notated_sequence = enz.value("notated_sequence", "");
        e.forward_sequence = enz.value("forward_sequence", "");
        e.reverse_complement = enz.value("reverse_complement", "");
        e.top_cuts = enz.value("top_cuts", std::vector<int>{});
        e.bottom_cuts = enz.value("bottom_cuts", std::vector<int>{});

        
        enzymes.push_back(e);
    }
    
    return enzymes;
}


/*
 * Calculates and returns the global nucleic acid index of the cut site of a given restriction enzyme
 * Takes as argument the find index of find_DNA_seq(RE_sequence, vector/gene_seq) and the local cut index of the restriction enzyme
 */
size_t determine_cut_index(size_t find_index, size_t local_RE_cut_index) {
    size_t global_cut_index;
    
    global_cut_index = find_index + local_RE_cut_index;
    
    return global_cut_index;
}
   

/*
 * Determines whether a string recognition_seq is found a SECOND time in a string main_sequence.
 * Assumes that the substr was found once at size_t index_found.
 */
bool found_again(const size_t& index_found, const std::string& recognition_seq, const std::string& main_seq) {
    size_t length_rec_seq = recognition_seq.size();
    // Rest of the main_seq after substr
    std::string new_main_seq = main_seq.substr(index_found + length_rec_seq);
    
    // If found again, return true
    if (nucleotide::find_DNA_seq(recognition_seq, new_main_seq) != std::string::npos) {
        return true;
    }
    
    return false;
}


/*
 * Find an available output filename given a base_name and extension.
 *
 */
std::string find_available_filename(const std::string& base_name, const std::string& extension) {
    int counter = 1;
    std::string filename = base_name + std::to_string(counter)+ extension;

    // Keep checking for existence and incrementing the counter until a unique name is found
    while (std::filesystem::exists(filename)) {
        counter++;
        filename = base_name + std::to_string(counter) + extension;
    }

    return filename;
}

/*
 * Writes a linear cut map of a vector to a .txt file. Shows where restriction enzymes cut the vector.
 */
void write_linear_cut_map(
    const std::string vector_name, 
    const std::string &vector_seq, 
    const std::vector<RestrictionEnzyme> &enzymes
) {
    constexpr size_t wrap_width = 80;
    const size_t vector_length = vector_seq.size();
    const size_t set_w_1 = 20;
    const size_t set_w_2 = 15;

    std::string file_name = find_available_filename("output/linear_map", ".txt");
    std::ofstream out(file_name);
    if (!out) {
        throw std::runtime_error("Failed to open " + file_name + "for writing");
    }
    
    /* Informational heading. */
    out << " - LINEAR MAP OF VECTOR CUT SITES -\n\n";
    out << "Vector file: " << vector_filename << "\n"
        << "Gene file: " << gene_filename << "\n\n";

    out << std::left
        << std::setw(set_w_1) << "Enzyme"
        << std::setw(set_w_2) << "Site"
        << "Cut position(s)\n";
    out << "--------------------------------------------------\n";

    // Sort enzymes alphabetically for readability
    std::vector<const RestrictionEnzyme *> sorted;
    for (const auto &e : enzymes) {
        sorted.push_back(&e);
    }

    std::sort(sorted.begin(), sorted.end(),
              [](const RestrictionEnzyme *a, const RestrictionEnzyme *b) {
                  return a->name < b->name;
              });

    for (const RestrictionEnzyme *enzyme : sorted) {
        out << std::left
            << std::setw(set_w_1) << enzyme->name
            << std::setw(set_w_2) << enzyme->notated_sequence;

        bool printed = false;

        if (enzyme->forward_RE_site_index != std::string::npos) {
            size_t cut_index = enzyme->forward_RE_site_index + enzyme->top_cuts[0]; // 1-based
            out << (cut_index);
            printed = true;
        }

        if (enzyme->forward_RE_site_index != std::string::npos && enzyme->revcomp_RE_site_index != enzyme->forward_RE_site_index) {
            if (printed) out << ", ";
            size_t cut_index = enzyme->revcomp_RE_site_index + (enzyme->forward_sequence.size() - enzyme->bottom_cuts[0]); // 1-based
            out << (cut_index);
        }

        out << '\n';
    }

    out << "\n\n" << vector_name << " (" << vector_length << " bp)" <<  "\n";
    
    
    /* Linear vector map. */

    // ---- Collect cuts ----
    struct Cut {
        // 0-based. Actual location of cut.
        // Ex.1: index = 3, vector_seq = "ATGCA"... means cut is between 3 and 4 -> ATGC/A
        size_t index;
        
        // Enzyme name.
        std::string name;
    };

    std::vector<Cut> cuts;

    for (const auto &e : enzymes) {
        if (e.forward_RE_site_index != std::string::npos && e.forward_RE_site_index < vector_length) {
            size_t cut_index = e.forward_RE_site_index + e.top_cuts[0] - 1;
            cuts.push_back({ cut_index, e.name });
        }
        if (e.revcomp_RE_site_index != std::string::npos && e.revcomp_RE_site_index < vector_length) {
            size_t cut_index = e.revcomp_RE_site_index + (e.forward_sequence.size() - e.bottom_cuts[0]) - 1;
            cuts.push_back({ cut_index, e.name });
        }
    }

    // Iterate over wrapped lines
    for (size_t line_start = 0; line_start < vector_length; line_start += wrap_width) {
        const size_t line_end = std::min(line_start + wrap_width, vector_length);
        const size_t line_len = line_end - line_start;

        /* Base pair numbering line */
        out << std::setw(10)<< (line_start + 1);
        for (size_t i = 10; i < line_len; i += 10) {
            out << std::setw(10) << (line_start + i + 1);
        }
        out << '\n';

        /* Sequence line */
        out << vector_seq.substr(line_start, line_len) << '\n';

        /* Annotation line */
        std::string annotation(line_len, ' ');

        for (const auto &cut : cuts) {
            if (cut.index < line_start || cut.index >= line_end)
                continue;

            const size_t caret_pos = cut.index - line_start;
            annotation[caret_pos] = '^';

            const size_t start = caret_pos + 1;
            if (start >= annotation.size())
                continue;

            // Measure available free space
            size_t available = 0;
            while (start + available < annotation.size() &&
                   annotation[start + available] == ' ') {
                ++available;
            }

            if (available == 0)
                continue;

            const std::string &name = cut.name;

            // Generate candidate labels (in priority order)
            std::vector<std::string> candidates;

            // 1. Full name
            candidates.push_back(" " + name);

            // 2. Strip common suffixes
            {
                std::string stripped = name;
                const std::vector<std::string> suffixes = { "-HF", "-v2", "-v3" };
                for (const auto &s : suffixes) {
                    if (stripped.size() > s.size() &&
                        stripped.compare(stripped.size() - s.size(), s.size(), s) == 0) {
                        stripped.erase(stripped.size() - s.size());
                        break;
                    }
                }
                if (stripped != name) {
                    candidates.push_back(" " + stripped);
                }
            }

            // 3–5. Progressive truncation
            if (name.size() >= 4)
                candidates.push_back(" " + name.substr(0, 4));
            if (name.size() >= 2)
                candidates.push_back(" " + name.substr(0, 2));
            if (!name.empty())
                candidates.push_back(" " + name.substr(0, 1));

            // Try candidates in order
            for (const auto &label : candidates) {
                if (label.size() <= available) {
                    for (size_t i = 0; i < label.size(); ++i) {
                        annotation[start + i] = label[i];
                    }
                    break;
                }
            }
        }

        out << annotation << "\n\n";
    }
    
    std::cout << "\nA linear map of cut sites has been output to: " << file_name.substr(7);
}

/*
 * Extracts name of a vector from a vector header file.
 */
std::string extract_vector_name(const std::string& fasta_header) {
    if (fasta_header.empty() || fasta_header[0] != '>') {
        throw std::invalid_argument("Invalid FASTA header");
    }

    // Remove leading '>'
    std::string header = fasta_header.substr(1);

    // Find first delimiter
    size_t pipe_pos = header.find('|');

    if (pipe_pos == std::string::npos) {
        return header;               // No metadata present
    }

    return header.substr(0, pipe_pos);
}

/*
 * Writes a circular vector/cut site map to an svg file.
 */
void write_svg_vector_map(
    const std::string &vector_name,
    const std::string &vector_seq,
    const std::vector<RestrictionEnzyme> &enzymes
) {
    namespace fs = std::filesystem;

    // ---- Output path ----
    std::filesystem::create_directories("output");
    std::string filename = find_available_filename("output/circular_map", ".svg");
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Failed to open " + filename + " for writing.");
    }

    // ---- Basic geometry ----
    const double W = 2100.0;
    const double H = 2100.0;
    const double cx = W / 2.0;
    const double cy = H / 2.0;

    const double r_backbone = 650.0;
    const double tick_len = 18.0;
    
    const size_t major_every = 8;     // 128 ticks / 16 = 8 major ticks
    const double major_tick_len = 28.0;
    const double major_tick_width = 4.0;
    
    const double tick_label_r = r_backbone - major_tick_len - 18.0; // inside ring

    const double cut_len_in = 14.0;
    const double cut_len_out = 22.0;

    const double label_r_even = r_backbone + 110.0;
    const double label_r_odd  = r_backbone + 170.0;

    const size_t L = vector_seq.size();
    if (L == 0) {
        throw std::runtime_error("vector_seq is empty; cannot draw vector map");
    }

    const double PI = 3.14159265358979323846;
    const double TWO_PI = 2.0 * PI;

    auto angle_for_bp = [&](size_t p) -> double {
        // 0-based bp index p in [0, L-1]
        double f = static_cast<double>(p) / static_cast<double>(L);
        double theta = TWO_PI * f;
        // rotate so bp0 is at 12 o'clock
        return theta - PI / 2.0;
    };

    auto xy_on_radius = [&](double theta, double r) -> std::pair<double, double> {
        return { cx + r * std::cos(theta), cy + r * std::sin(theta) };
    };

    // XML escaping for &,<,>
    auto escape_xml = [](const std::string &s) -> std::string {
        std::string r;
        r.reserve(s.size());
        for (char ch : s) {
            switch (ch) {
                case '&': r += "&amp;"; break;
                case '<': r += "&lt;"; break;
                case '>': r += "&gt;"; break;
                case '"': r += "&quot;"; break;
                case '\'': r += "&apos;"; break;
                default: r += ch; break;
            }
        }
        return r;
    };

    // ---- Collect cuts ----
    struct Cut {
        // 0-based. Actual location of cut. 
        // Ex.1: index = 3, vector_seq = "ATGCA"... means cut is between 3 and 4 -> ATGC/A
        size_t index;
        
        // Enzyme name.
        std::string name;
    };

    std::vector<Cut> cuts;
    cuts.reserve(enzymes.size() * 2);

    for (const auto &e : enzymes) {
        if (e.forward_RE_site_index != std::string::npos && e.forward_RE_site_index < L) {
            size_t cut_index = e.forward_RE_site_index + e.top_cuts[0] - 1;
            cuts.push_back({ cut_index, e.name });
        }
        if (e.revcomp_RE_site_index != std::string::npos && e.revcomp_RE_site_index < L) {
            size_t cut_index = e.revcomp_RE_site_index + (e.forward_sequence.size() - e.bottom_cuts[0]) - 1;
            cuts.push_back({ cut_index, e.name });
        }
    }

    std::sort(cuts.begin(), cuts.end(), [](const Cut &a, const Cut &b) {
        if (a.index != b.index) return a.index < b.index;
        return a.name < b.name;
    });

    // Dedupe exact duplicates (same index + name)
    cuts.erase(std::unique(cuts.begin(), cuts.end(), [](const Cut &a, const Cut &b) {
        return a.index == b.index && a.name == b.name;
    }), cuts.end());

    // ---- Identify dense regions and compute per-cut label angle/radius overrides ----
    struct Override {
        bool on = false;
        double theta = 0.0; // label angle (radians)
        double r = 0.0;     // label radius
    };
    std::vector<Override> ov(cuts.size());

    if (cuts.size() >= 4) {
        // Heuristics:
        // "Dense" = runs of consecutive cuts where gap <= avg_spacing/3 (min 12 bp), and cluster size >= 4.
        const double avg_spacing = static_cast<double>(L) / static_cast<double>(cuts.size());
        const double gap_threshold_bp = std::max(100.0, avg_spacing / 3.0);
        const size_t min_cluster_size = 4;

        // Leader radius range for dense region (midpoint longest)
        const double dense_r_max = r_backbone + 320.0; // longest (middle)
        const double dense_r_min = r_backbone + 120.0; // shortest (edges)

        // Max angular fan-out at cluster edges
        const double max_fan_deg = 9.0;
        const double max_fan_rad = (PI / 180.0) * max_fan_deg;

        struct Cluster {
            size_t start_k; // inclusive
            size_t end_k;   // inclusive
            bool wraps = false;
        };

        std::vector<Cluster> clusters;
        clusters.reserve(cuts.size() / 2);

        size_t start = 0;
        for (size_t k = 1; k < cuts.size(); ++k) {
            double gap = static_cast<double>(cuts[k].index) - static_cast<double>(cuts[k - 1].index);
            if (gap <= gap_threshold_bp) continue;

            if (k - start >= min_cluster_size) {
                clusters.push_back({ start, k - 1, false });
            }
            start = k;
        }
        if (cuts.size() - start >= min_cluster_size) {
            clusters.push_back({ start, cuts.size() - 1, false });
        }

        // Merge wrap-around cluster if needed
        if (!clusters.empty()) {
            double wrap_gap = (static_cast<double>(L) - static_cast<double>(cuts.back().index)) +
                              static_cast<double>(cuts.front().index);

            const bool last_is_tail_cluster = (clusters.back().end_k == cuts.size() - 1);
            const bool first_is_head_cluster = (clusters.front().start_k == 0);

            if (last_is_tail_cluster && first_is_head_cluster && wrap_gap <= gap_threshold_bp) {
                Cluster merged;
                merged.start_k = clusters.back().start_k;
                merged.end_k = clusters.front().end_k;
                merged.wraps = true;

                clusters.erase(clusters.begin());
                clusters.pop_back();
                clusters.push_back(merged);
            }
        }

        // Apply overrides for each cluster
        for (const auto &cl : clusters) {
            std::vector<size_t> ks; // indices into cuts for this cluster (in cluster order)
            ks.reserve(cuts.size());

            if (!cl.wraps) {
                for (size_t k = cl.start_k; k <= cl.end_k; ++k) ks.push_back(k);
            } else {
                for (size_t k = cl.start_k; k < cuts.size(); ++k) ks.push_back(k);
                for (size_t k = 0; k <= cl.end_k; ++k) ks.push_back(k);
            }

            if (ks.size() < min_cluster_size) continue;

            // "Middle": cut closest to cluster midpoint in terms of cut-count (index in ks)
            // (If even size, this picks the lower middle; that's fine.)
            const size_t mid_i = ks.size() / 2;
            const size_t mid_k = ks[mid_i];

            const double mid_theta = angle_for_bp(cuts[mid_k].index);

            // Distance-from-middle is based on cut-site count, not bp distance
            const double max_dist = static_cast<double>(std::max(mid_i, (ks.size() - 1) - mid_i));
            if (max_dist <= 0.0) {
                // cluster size 1 (shouldn't happen due to min_cluster_size)
                ov[mid_k] = { true, mid_theta, dense_r_max };
                continue;
            }

            for (size_t i = 0; i < ks.size(); ++i) {
                const size_t k = ks[i];
                const int d = static_cast<int>(i) - static_cast<int>(mid_i); // negative = before middle
                const double t = std::min(1.0, std::fabs(static_cast<double>(d)) / max_dist);

                // Longest at middle, shorter toward edges
                const double r_here = dense_r_max - (dense_r_max - dense_r_min) * t;

                // Lean away from middle: sign depends on side, magnitude grows with |d|
                const double sign = (d < 0) ? -1.0 : (d > 0 ? 1.0 : 0.0);
                const double theta_here = angle_for_bp(cuts[k].index) + sign * max_fan_rad * t;

                ov[k].on = true;
                ov[k].theta = theta_here;
                ov[k].r = r_here;
            }

            // Enforce exact middle behavior
            ov[mid_k].on = true;
            ov[mid_k].theta = mid_theta;
            ov[mid_k].r = dense_r_max;
        }
    }

    // ---- SVG header ----
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
        << "style=\"width:100vw; height:100vh; display:block;\" "
        << "viewBox=\"0 0 " << static_cast<int>(W) << " " << static_cast<int>(H) << "\" "
        << "preserveAspectRatio=\"xMidYMid meet\">\n";

    // Styling
    out << "  <style>\n"
        << "    .backbone { fill: none; stroke: #000; stroke-width: 6; }\n"
        << "    .tick { stroke: #000; stroke-width: 2; }\n"
        << "    .cut { stroke: #000; stroke-width: 3; }\n"
        << "    .leader { stroke: #000; stroke-width: 2; }\n"
        << "    .label { font-family: sans-serif; font-size: 22px; }\n"
        << "    .center1 { font-family: sans-serif; font-size: 34px; font-weight: 600; }\n"
        << "    .center2 { font-family: sans-serif; font-size: 26px; }\n"
        << "    .tick_label { font-family: sans-serif; font-size: 18px; }\n"
        << "  </style>\n";

    // White background
    out << "  <rect x=\"0\" y=\"0\" width=\"" << static_cast<int>(W)
        << "\" height=\"" << static_cast<int>(H) << "\" fill=\"white\" />\n";

    // ---- Backbone ----
    out << "  <circle class=\"backbone\" cx=\"" << cx << "\" cy=\"" << cy
        << "\" r=\"" << r_backbone << "\" />\n";

    // ---- Tick marks: exactly 128 ----
    const size_t tick_count = 128;
    for (size_t i = 0; i < tick_count; ++i) {
        double theta = (2.0 * PI * static_cast<double>(i) / static_cast<double>(tick_count)) - PI / 2.0;

        const bool is_major = (i % major_every == 0);
        const double len = is_major ? major_tick_len : tick_len;
        const double w = is_major ? major_tick_width : 2.0;

        auto [x1, y1] = xy_on_radius(theta, r_backbone - len);
        auto [x2, y2] = xy_on_radius(theta, r_backbone);

        out << "  <line class=\"tick\" x1=\"" << x1 << "\" y1=\"" << y1
            << "\" x2=\"" << x2 << "\" y2=\"" << y2
            << "\" stroke-width=\"" << w << "\" />\n";

        if (is_major) {
            // map major tick index to bp position (1-based for display)
            const size_t bp = (i * L) / tick_count;   // 0..L-1
            const size_t bp_display = bp + 1;

            auto [tx, ty] = xy_on_radius(theta, tick_label_r);

            // anchor radially
            const std::string anchor = "middle";

            const double radial_pad = 2.0;

            double lx = tx + radial_pad * std::cos(theta);
            double ly = ty + radial_pad * std::sin(theta);
            
            // If we're near 3 o'clock or 9 o'clock, pull the label inward to avoid tick overlap.
            const double near_side = 0.20;          // smaller = stricter "only the two sides"
            const double side_inset = 10.0;         // px inward; tune this

            if (std::fabs(std::sin(theta)) < near_side) {
                lx -= side_inset * std::cos(theta);
                ly -= side_inset * std::sin(theta);
            }
            
            out << "  <text class=\"tick_label\" x=\"" << lx << "\" y=\"" << ly
                << "\" text-anchor=\"" << anchor << "\" dominant-baseline=\"middle\">"
                << bp_display << "</text>\n";
        }

    }

    // ---- Cut marks ----
    for (const auto &c : cuts) {
        double theta = angle_for_bp(c.index);

        auto [x1, y1] = xy_on_radius(theta, r_backbone - cut_len_in);
        auto [x2, y2] = xy_on_radius(theta, r_backbone + cut_len_out);

        out << "  <line class=\"cut\" x1=\"" << x1 << "\" y1=\"" << y1
            << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" />\n";
    }

    // ---- Labels (alternating radius + dense overrides) ----
    for (size_t k = 0; k < cuts.size(); ++k) {
        const auto &c = cuts[k];

        // True cut angle for leader start (always correct)
        const double theta_cut = angle_for_bp(c.index);

        // Defaults
        double theta_label = theta_cut;
        double r_label = (k % 2 == 0) ? label_r_even : label_r_odd;

        // Dense override
        if (k < ov.size() && ov[k].on) {
            theta_label = ov[k].theta;
            r_label = ov[k].r;
        }

        auto [lx1, ly1] = xy_on_radius(theta_cut, r_backbone + cut_len_out + 6.0);
        auto [lx2, ly2] = xy_on_radius(theta_label, r_label - 10.0);

        out << "  <line class=\"leader\" x1=\"" << lx1 << "\" y1=\"" << ly1
            << "\" x2=\"" << lx2 << "\" y2=\"" << ly2 << "\" />\n";

        bool right_side = std::cos(theta_label) >= 0.0;
        std::string anchor = right_side ? "start" : "end";
        double pad = 8.0 * (right_side ? 1.0 : -1.0);

        std::string label = escape_xml(c.name + " (" + std::to_string(c.index + 1) + ")");

        out << "  <text class=\"label\" x=\"" << (lx2 + pad) << "\" y=\"" << ly2
            << "\" text-anchor=\"" << anchor << "\" dominant-baseline=\"middle\">"
            << label << "</text>\n";
    }

    // ---- Center title ----
    out << "  <text class=\"center1\" x=\"" << cx << "\" y=\"" << (cy - 18)
        << "\" text-anchor=\"middle\" dominant-baseline=\"middle\">"
        << escape_xml(vector_name) << "</text>\n";

    out << "  <text class=\"center2\" x=\"" << cx << "\" y=\"" << (cy + 22)
        << "\" text-anchor=\"middle\" dominant-baseline=\"middle\">"
        << L << " bp</text>\n";

    out << "</svg>\n";

    std::cout << "\nA circular map of cut sites has been output to: " << filename.substr(7);
}



/*
 * Queries user until the user gives a valid input. Returns valid user input.
 * If any input is valid, param valid_inputs should be input as empty vector
 */
std::string query_user (
    const std::string& query,
    const std::vector<std::string>& valid_inputs,
    const std::string& invalid_input_message
) {
    std::string user_input;
    
    if (valid_inputs.size() == 0) {
        std::cout << "\n" << query << ": ";
        std::getline(std::cin, user_input);
        return user_input;
    }
    
    bool good_input = false;
    while (!good_input) {
        std::cout << "\n" << query << ": ";
        std::getline(std::cin, user_input);
        for (std::string vi: valid_inputs) {
            if (vi.find("#length:") != std::string::npos) {
                if (std::to_string(user_input.size()) == std::string (1, vi[vi.find(":")+1])) { // very dirty
                    good_input = true;
                }
            } else if (user_input == vi) {
                good_input = true;
            }
            
        }
        if (!good_input) {
            std::cout << "\n" << invalid_input_message << "\n";
        }
    }
    return user_input;
}


/*
 * Generates DNA sequences that code for given protein tags.
 * tag_id has 2 uses: it either contains the name of the tag (ex. "His6") OR the tag DNA sequence (if it does not
 *  contain any built-in tag names).
 */
std::string generate_tag_seq(std::string& tag_id) {
    std::string tag_seq = "";
    if (tag_id == "His6") {
        // HHHHHH
        std::string his2_seq = "CATCAC";
        tag_seq = repeat(his2_seq, 3);   //CATCACCATCACCATCAC
        return tag_seq;
        
    } else if (tag_id == "FLAG") {
        // DYKDDDDK
        tag_seq = "GATTATAAAGATGATGATGATAAATAA";
        return tag_seq;
        
    } else if (tag_id == "HA") {
        // YPYDVPDYA
        tag_seq = "TATCCGTATGATGTGCCGGATTATGCGTAA";
        return tag_seq;
        
    } else if (tag_id == "Myc") {
        // EQKLISEEDL
        tag_seq = "GAACAGAAACTGATTAGCGAAGAAGATCTGTAA";
        return tag_seq;
        
    } else { // tag_id IS the DNA sequence
        return tag_id;
    }
}

/*
 * Generates DNA sequences that code for given protease recognition sites.
 * protease_id has 2 uses: it either contains the name of the protease (ex. "TEV") OR the protease DNA sequence (if it
 *  does not contain any built-in protease names).
 */
std::string generate_protease_seq(std::string& protease_id) {
    std::string protease_seq = "";
    if (protease_id == "TEV") {
        // ENLYFQS
        protease_seq = "GAAAACCTGTATTTTCAGAGCTAA";
        return protease_seq;

    } else if (protease_id == "PreScission") {
        // LEVLFQGP
        protease_seq = "CTGGAAGTGCTGTTTCAGGGCCCGTAA";
        return protease_seq;
        
    } else if (protease_id == "FXa") {
        // IEDGR
        protease_seq = "ATTGAAGATGGCCGCTAA";
        return protease_seq;
        
    } else if (protease_id == "Enterokinase") {
        // DDDDK
        protease_seq = "GATGATGATGATAAATAA";
        return protease_seq;
        
    } else { // protease_id IS the DNA recognition sequence
        return protease_id;
    }
}


/*
 * Returns str repeated num_times times.
 */
std::string repeat(std::string str, size_t num_times) {
    std::string out = "";
    if (num_times > std::string::npos-50) return "Error: too many repeats";
    for (size_t i = 0; i < num_times; i++) {
        out += str;
    }
    return out;
}


/* ---------------------------------- END FUNCTION DEFINITIONS --------------------------------- */
/* --------------------------------------------------------------------------------------------- */

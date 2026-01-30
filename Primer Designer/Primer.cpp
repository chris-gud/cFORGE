// Primer.cpp
#include "Primer.hpp"
#include "nucleotide_utils.hpp"

/*
 * Primer constructor (std::string sequence)
 */
Primer::Primer (std::string sequence) {
    // Assign values to seq, length, and wallace_melting_temp
    seq = sequence;
    length = seq.size();
            
    // Calculate number of each amino acid and assign to corresponding variable
    for (char c : seq) {
        if (c == 'A')      num_A++;
        else if (c == 'T') num_T++;
        else if (c == 'G') num_G++;
        else if (c == 'C') num_C++;
    }
    
    // Calculate Wallace melting temperature
    wallace_melting_temp = 2*num_A + 2*num_T + 4*num_C + 4*num_G;
    
    // Calculate SantaLucia melting temperature
    santalucia_melting_temp = calculate_tm_santalucia_2004(seq);
    
    // Calculate G and C content and assign to GC_content
    GC_content = (static_cast<float>(num_G) + static_cast<float>(num_C)) / (static_cast<float>(length));
    
    // Assess GC clamp within 5 nucleotides of 3' end
    std::string last_five = seq.substr(length-5);
    if (last_five.find('C') == std::string::npos && last_five.find('G') == std::string::npos)
        GC_clamp = false;
    else GC_clamp = true;
    
    // Assess whether primer undergoes homo dimerization
    homo_dimerizes = dimerizes_with(seq);
    
    // Determine the longest homopolymer run in seq
    max_length_homopolymer_run = 1;
    size_t length_run = 1;
    for (size_t i = 0; i < length - 1; i++) {
        if (seq.at(i) == seq.at(i+1)) {
            length_run++;
        } else {
            if (length_run > max_length_homopolymer_run) {
                max_length_homopolymer_run = length_run;
            }
            length_run = 1;
        }
    }
    
    // Assess whether primer forms deleterious hairpins
    hairpin = forms_hairpin();
}


/* Useful public Primer functions */

/*
 * Returns true if Primer this and Primer other are LIKELY to dimerize.
 * Good for reuse to determine homodimerization OR heterodimerization.
 * Utilizes the dimerize() function in nucleotide_utils.
 * returns true if dimerize(this->seq, other_seq) returns a double >= 0.75
 */
bool Primer::dimerizes_with(const std::string& other_seq) {
    double dimerization_prob = nucleotide::dimerize(this->seq, other_seq);
    return dimerization_prob >= 0.75;
}


/*
 * bool is_good() returns true if this Primer satisfies the primer design rules as notated below:
 *
 * Primer rules:
 *   1.) 18-25 bases long ideally. PCR mutagenesis allows up to 45 n.t.
 *   2.) 40-70% GC content
 *   3.) Melting temperature between 55-65°C
 *   4.) GC clamp at the 3' end
 *   5.) Avoid primer homodimerization, homopolymer runs, hairpins, off-target binding, heterodimerization
 *
 */
bool Primer::is_good() {
    /* ASSESS RULE #1: Sequence is 18-45 bases long */
    if (! (length >= 18 && length <= 45)) return false;
    
    /* ASSESS RULE #2: 40-70% GC content */
    if (! (GC_content >= 0.40 && GC_content <= 0.70)) return false;
    
    /* ASSESS RULE #3: SantaLucia melting temperature between 55-65°C */
    if (! (santalucia_melting_temp >= 55 && santalucia_melting_temp <= 65)) return false;
    
    /* ASSESS RULE #4: GC clamp at the 3' end */
    if (!GC_clamp) return false;
    
    /* ASSESS RULE #5A: Avoid primer homodimerization */
    if (homo_dimerizes) return false;
    
    /* ASSESS RULE #5B: Avoid homopolymer runs longer than 5 n.t. */
    size_t acceptable_length_homopolymer_run = 5;
    if (max_length_homopolymer_run > acceptable_length_homopolymer_run) return false;
    
    /* ASSESS RULE #5C: Avoid deleterious hairpins */
    if (hairpin) return false;
    
    /* RULES #5D & #5E CANNOT be assessed within this function. They will be assessed via public method in main. */
    
    return true;
}


/*
 * bool binds_off_target(std::string &gene_seq) returns true if the reverse complement of the
 *  primer sequence is found exactly in gene_seq OR the reverse complement of the gene_seq MORE THAN ONCE.
 * For use within main.
 * SHOULD BE improved to use a nearest neighbor model rather than naively comparing sequences.
 *  I will do that another time.
 */
bool Primer::binds_off_target(const std::string &gene_seq) {
    // Total number of primer binding sites within the forward AND reverse strands of the gene
    // NOTE: not ACTUALLY - since we only care about the first 2 on either strand, the maximum
    //  of total_num_binds is 4
    int total_num_binds = 0;
    // Generate reverse complement of primer (i.e. binding site of primer)
    std::string rev_comp_primer = nucleotide::revcomp_of(seq);
    
    /*
     Find number of binding sites between primer and FORWARD gene strand
     */
    // Find first location primer binds to on forward strand
    size_t index_first = gene_seq.find(rev_comp_primer);
    // Make sure new substring isn't out of bounds (in the case the only binding site is
    //  the proper reverse primer binding site)
    if (index_first != std::string::npos) {
        // There is at least one binding site. Increment
        total_num_binds++;
        if (index_first + rev_comp_primer.size() < gene_seq.size()) {
            // Generate substring of gene sequence minus the first.
            std::string gene_seq_substr = gene_seq.substr(index_first + rev_comp_primer.size());
            // If found again, return true
            if (gene_seq_substr.find(rev_comp_primer) != std::string::npos)
                // Another binding site was found; increment
                total_num_binds++;
        }
    }
    
    /*
     Find number of binding sites between primer and REVERSE gene strand
     */
    // Generate reverse complement of gene sequence (i.e. the reverse strand)
    std::string rev_comp_gene_seq = nucleotide::revcomp_of(gene_seq);
    // Find first location primer binds to on reverse strand
    size_t index_first_rev = rev_comp_gene_seq.find(rev_comp_primer);
    if (index_first_rev != std::string::npos) {
        // There is at least one binding site. Increment
        total_num_binds++;
        if (index_first_rev + rev_comp_primer.size() < rev_comp_gene_seq.size()) {
            // Generate substring of gene sequence minus the first.
            std::string gene_seq_substr_rev = rev_comp_gene_seq.substr(index_first_rev + rev_comp_primer.size());
            // If found again, return true
            if (gene_seq_substr_rev.find(rev_comp_primer) != std::string::npos)
                // Another binding site was found; increment
                total_num_binds++;
        }
    }
    
    // If the primer binds to the gene at >1 location, it binds off target
    if (total_num_binds > 1) {        
        return true;
    }
    
    return false;
}


/*
 * to_string() returns the Primer information in std::string format
 */
std::string Primer::to_string() {
    std::string string_seq = "5'- " + seq + " -3'";
    std::string string_length = std::to_string(length) + " n.t.";
    std::string string_wallace_melting_temp = std::to_string(wallace_melting_temp) + "°C";
    std::string string_sl_melting_temp = std::to_string(santalucia_melting_temp).substr(0, 4) + "°C";
    std::string string_GC_content = std::to_string(GC_content*100.).substr(0, 4) + "%";
    std::string string_GC_clamp = GC_clamp ? "YES" : "NO";
    std::string string_homo_dimerizes = homo_dimerizes ? "YES" : "NO";
    std::string string_max_length_homopolymer_run = std::to_string(max_length_homopolymer_run) + " n.t.";
    std::string string_hairpin = hairpin ? "YES" : "NO";
   
    std::string out = "";
    out += "----------------------------------------------------------------------\n";
    out += "        Primer Sequence: " + string_seq + "\n";
    out += "                 Length: " + string_length + "\n";
    out += "   Wallace Melting Temp: " + string_wallace_melting_temp + "\n";
    out += "SantaLucia Melting Temp: " + string_sl_melting_temp + "\n";
    out += "             GC Content: " + string_GC_content + + "\n";
    out += "               GC Clamp: " + string_GC_clamp + "\n";
    out += "       Homodimerization: " + string_homo_dimerizes + "\n";
    out += " Homopolymer Run Length: " + string_max_length_homopolymer_run + "\n";
    out += "       Forms Hairpin(s): " + string_hairpin + "\n";
    out += "----------------------------------------------------------------------\n";

    return out;
}


/*
 * bool forms_hairpin returns true if Primer sequence is likely to form a deleterious secondary structure.
 * Utilizes the hairpin() function in nucleotide_utils.
 * returns true if hairpin(this->seq) returns a double >= 0.65
 */
bool Primer::forms_hairpin() {
    double hairpin_prob = nucleotide::hairpin(this->seq);
    return hairpin_prob >= 0.65;
}


/*
 * Calculates the melting temperature of the Primer using SantaLucia 2004.
 */
double Primer::calculate_tm_santalucia_2004(
    const std::string& primer,
    double primer_conc,    // 500 nM default
    double na_conc         // 50 mM Na+ default
) {
    if (primer.size() < 2)
        throw std::invalid_argument("Primer too short");
    
    // Nearest-neighbor parameters (SantaLucia 2004)
    // ΔH (kcal/mol), ΔS (cal/mol·K)
    static const std::unordered_map<std::string, std::pair<double, double>> nn = {
        {"AA", {-7.9, -22.2}}, {"TT", {-7.9, -22.2}},
        {"AT", {-7.2, -20.4}}, {"TA", {-7.2, -21.3}},
        {"CA", {-8.5, -22.7}}, {"TG", {-8.5, -22.7}},
        {"GT", {-8.4, -22.4}}, {"AC", {-8.4, -22.4}},
        {"CT", {-7.8, -21.0}}, {"AG", {-7.8, -21.0}},
        {"GA", {-8.2, -22.2}}, {"TC", {-8.2, -22.2}},
        {"CG", {-10.6, -27.2}},
        {"GC", {-9.8, -24.4}},
        {"GG", {-8.0, -19.9}}, {"CC", {-8.0, -19.9}}
    };
    
    double dH = 0.0;
    double dS = 0.0;
    
    auto base = [&](char c) {
        c = std::toupper(c);
        if (c != 'A' && c != 'T' && c != 'G' && c != 'C')
            throw std::invalid_argument("Invalid nucleotide in primer");
        return c;
    };
    
    // Sum nearest neighbors
    for (size_t i = 0; i + 1 < primer.size(); i++) {
        std::string pair;
        pair += base(primer[i]);
        pair += base(primer[i + 1]);
        
        auto it = nn.find(pair);
        if (it == nn.end())
            throw std::runtime_error("NN lookup failure");
        
        dH += it->second.first;
        dS += it->second.second;
    }
    
    // Initiation parameters
    dH += 0.2;
    dS += -5.7;
    
    // Terminal AT penalty
    char first = base(primer.front());
    char last  = base(primer.back());
    
    if (first == 'A' || first == 'T') {dH += 2.2; dS += 6.9;}
    if (last  == 'A' || last  == 'T') {dH += 2.2; dS += 6.9;}

    // Gas constant
    const double R = 1.987; // cal/mol·K

    // Melting temperature (Kelvin)
    double tm_kelvin =
        (dH * 1000.0) /
        (dS + R * std::log(primer_conc / 4.0));

    // Convert to Celsius and apply salt correction
    double tm_celsius =
        tm_kelvin - 273.15 + 16.6 * std::log10(na_conc);

    return tm_celsius;
}

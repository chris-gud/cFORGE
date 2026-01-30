// RestrictionEnzyme.cpp
#include "RestrictionEnzyme.hpp"
#include "nucleotide_utils.hpp"

std::string RestrictionEnzyme::to_string() {
    std::string out = "";
    out += "--------------------------------------------------------------\n";
    out += "Restriction Enzyme Name: " + name + "  |  Notated Sequence: " + notated_sequence;
    out += "\nForward Sequence: " + forward_sequence + "  |  Reverse Complement: " + reverse_complement;
    
    out += "\nTop Cuts: ";
    for (int cut : top_cuts)
        out += std::to_string(cut) + " ";
    out += "  |  Bottom Cuts: ";
    for (int cut : bottom_cuts)
        out += std::to_string(cut) + " ";
    
    std::string string_forward_cuts_vector = forward_cuts_vector ? "YES" : "NO";
    std::string string_revcomp_cuts_vector = revcomp_cuts_vector ? "YES" : "NO";
    std::string string_forward_cuts_gene = forward_cuts_gene ? "YES" : "NO";
    std::string string_revcomp_cuts_gene = revcomp_cuts_gene ? "YES" : "NO";

    
    out += "\nForward cuts vector: " + string_forward_cuts_vector + "  |  Revcomp cuts vector: " + string_revcomp_cuts_vector;
    out += "\nForward cuts gene: " + string_forward_cuts_gene + "  |  Revcomp cuts gene: " + string_revcomp_cuts_gene;
    out += "\nForward sequence index: " + std::to_string(forward_RE_site_index) + "  |  Revcomp sequence index: " + std::to_string(revcomp_RE_site_index);
    out += "\n";

    return out;
}

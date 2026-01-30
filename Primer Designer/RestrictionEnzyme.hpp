// RestrictionEnzyme.hpp
#ifndef RestrictionEnzyme_h
#define RestrictionEnzyme_h

#include <string>

struct RestrictionEnzyme {
    /* Data from JSON */
    std::string name;
    std::string notated_sequence; // Ex. (10/15)ACNNNNGTAYC(12/7) Ex. A/CCGGT
    std::string forward_sequence; // Forward sequence of the recognition site (5' -> 3')
    std::string reverse_complement; // Reverse complement of the recognition site (5' -> 3')
    std::vector<int> top_cuts;        // Local cut indices 5' -> 3' of top strand (can be negative)
    std::vector<int> bottom_cuts;     // Local cut indices 5' -> 3' of top strand (can be negative)
    /* End data from JSON */
    
    bool forward_cuts_vector = false;
    bool revcomp_cuts_vector = false;

    bool forward_cuts_gene = true;   // Default is true (since we are looking for falses)
    bool revcomp_cuts_gene = true;
    
    size_t forward_RE_site_index = std::string::npos;   // location of sequence in vector (0-based)
    size_t revcomp_RE_site_index = std::string::npos;   // location of revcomp sequence in vector (0-based)
        
    std::string to_string();
};



#endif /* RestrictionEnzyme_h */

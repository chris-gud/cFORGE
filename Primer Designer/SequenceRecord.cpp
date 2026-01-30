#include "SequenceRecord.hpp"

// ---------------- Feature helpers ----------------

std::optional<std::string> SequenceRecord::Feature::qualifier_first(const std::string& qkey) const {
    auto it = qualifiers.find(qkey);
    if (it == qualifiers.end() || it->second.empty()) return std::nullopt;
    return it->second.front();
}

void SequenceRecord::Feature::add_qualifier(const std::string& qkey, const std::string& value) {
    qualifiers[qkey].push_back(value);
}

// ---------------- SequenceRecord core ----------------

SequenceRecord::SequenceRecord(std::string name, std::string sequence)
    : name_(std::move(name)) {
    set_sequence(std::move(sequence));
}

void SequenceRecord::set_sequence(std::string seq) {
    sequence_ = normalize_sequence(seq);
}

void SequenceRecord::set_raw_field(std::string key, std::string value) {
    raw_fields_[std::move(key)] = std::move(value);
}

void SequenceRecord::add_feature(Feature f) {
    features_.push_back(std::move(f));
}

std::vector<const SequenceRecord::Feature*> SequenceRecord::find_features_by_key(const std::string& key) const {
    std::vector<const Feature*> out;
    for (const auto& f : features_) {
        if (f.key == key) out.push_back(&f);
    }
    return out;
}

bool SequenceRecord::is_valid_dna() const {
    if (sequence_.empty()) return false;
    for (char c : sequence_) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') return false;
    }
    return true;
}

// ---------------- Loading by path ----------------

SequenceRecord SequenceRecord::load_fasta_file(const std::string& path) {
    std::ifstream fin(path);
    if (!fin) throw std::runtime_error("SequenceRecord::load_fasta_file: could not open " + path);
    return load_fasta(fin);
}

SequenceRecord SequenceRecord::load_genbank_file(const std::string& path) {
    std::ifstream fin(path);
    if (!fin) throw std::runtime_error("SequenceRecord::load_genbank_file: could not open " + path);
    return load_genbank(fin);
}

static bool is_fasta_stream(std::istream& in) {
    // Peek first non-whitespace character
    std::streampos pos = in.tellg();
    in.clear();
    char ch = 0;
    while (in.get(ch)) {
        if (!std::isspace((unsigned char)ch)) break;
    }
    in.clear();
    in.seekg(pos);
    return ch == '>';
}

SequenceRecord SequenceRecord::load_from_file(const std::string& path) {
    std::ifstream fin(path);
    if (!fin) throw std::runtime_error("SequenceRecord::load_from_file: could not open " + path);

    if (is_fasta_stream(fin)) return load_fasta(fin);
    return load_genbank(fin);
}

// ---------------- Utility helpers ----------------

std::string SequenceRecord::trim(const std::string& s) {
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) ++a;
    while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}

std::string SequenceRecord::normalize_sequence(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        unsigned char uc = (unsigned char)c;
        if (std::isspace(uc)) continue;
        if (std::isdigit(uc)) continue; // useful for ORIGIN blocks with numbers
        if (!std::isalpha(uc)) continue;

        char up = (char)std::toupper(uc);
        if (up == 'U') up = 'T'; // treat RNA as DNA for primer design
        out.push_back(up);
    }
    return out;
}

static bool starts_with(const std::string& s, const std::string& pfx) {
    return s.size() >= pfx.size() && s.compare(0, pfx.size(), pfx) == 0;
}

static std::string slice_cols(const std::string& line, std::size_t start, std::size_t len = std::string::npos) {
    if (start >= line.size()) return "";
    if (len == std::string::npos) return line.substr(start);
    return line.substr(start, std::min(len, line.size() - start));
}

// ---------------- FASTA parsing ----------------

SequenceRecord SequenceRecord::load_fasta(std::istream& in) {
    SequenceRecord v;

    std::string line;
    bool saw_header = false;
    std::ostringstream seq;

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (saw_header) {
                // For simplicity: single-record FASTA. Stop at next header.
                break;
            }
            saw_header = true;
            std::string header = trim(line.substr(1));
            // Store whole header as name; you can split at first space if you want.
            v.name_ = header;
            continue;
        }

        // Sequence line
        seq << line;
    }

    v.set_sequence(seq.str());

    // If no name, give a fallback
    if (v.name_.empty()) v.name_ = "FASTA_record";

    // Topology unknown for FASTA
    v.topology_ = Topology::Unknown;

    return v;
}

// ---------------- GenBank parsing ----------------
//
// This is a pragmatic GenBank flatfile parser:
// - Captures common header fields
// - Parses FEATURES entries into Feature {key, location, qualifiers}
// - Captures ORIGIN sequence
// - Stores unknown header keys into raw_fields_
//
// Location expressions are stored as raw strings (no join/complement parsing here).

SequenceRecord SequenceRecord::load_genbank(std::istream& in) {
    enum class Section { Header, Features, Origin, Done };
    Section section = Section::Header;

    SequenceRecord v;

    std::string line;
    std::string current_header_key;

    // Feature state
    Feature* current_feature = nullptr;

    // Qualifier multiline state
    bool in_qual_value = false;
    std::string q_key;
    std::string q_value;
    bool q_quoted = false;

    auto commit_multiline_qual = [&]() {
        if (!q_key.empty() && current_feature) {
            current_feature->add_qualifier(q_key, q_value);
        }
        in_qual_value = false;
        q_key.clear();
        q_value.clear();
        q_quoted = false;
    };

    auto append_opt = [&](std::optional<std::string>& dst, const std::string& more) {
        std::string t = trim(more);
        if (t.empty()) return;
        if (!dst.has_value()) dst = t;
        else dst = *dst + " " + t;
    };

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();

        if (starts_with(line, "//")) {
            if (in_qual_value) commit_multiline_qual();
            section = Section::Done;
            break;
        }

        if (section == Section::Header) {
            if (starts_with(line, "FEATURES")) {
                section = Section::Features;
                current_header_key.clear();
                continue;
            }
            if (starts_with(line, "ORIGIN")) {
                section = Section::Origin;
                current_header_key.clear();
                continue;
            }

            // GenBank header: key typically in cols 1-12; value starts at col 13
            std::string key = trim(slice_cols(line, 0, 12));
            std::string val = slice_cols(line, 12);

            if (!key.empty()) {
                current_header_key = key;

                if (key == "LOCUS") {
                    v.locus_ = trim(val);
                    // Try to infer topology from LOCUS line text
                    std::string locus_line = v.locus_.value_or("");
                    std::string up = locus_line;
                    std::transform(up.begin(), up.end(), up.begin(),
                                   [](unsigned char c){ return (char)std::toupper(c); });
                    if (up.find("CIRCULAR") != std::string::npos) v.topology_ = Topology::Circular;
                    else if (up.find("LINEAR") != std::string::npos) v.topology_ = Topology::Linear;
                } else if (key == "DEFINITION") {
                    append_opt(v.definition_, val);
                } else if (key == "ACCESSION") {
                    v.accession_ = trim(val);
                } else if (key == "VERSION") {
                    v.version_ = trim(val);
                } else if (key == "SOURCE") {
                    append_opt(v.source_, val);
                } else if (key == "ORGANISM") {
                    append_opt(v.organism_, val);
                } else {
                    // Keep everything else for forward compatibility
                    // For multi-line values, we'll append on continuation lines below.
                    v.raw_fields_[key] = trim(val);
                }
            } else if (!current_header_key.empty()) {
                // Continuation line
                if (current_header_key == "DEFINITION") append_opt(v.definition_, val);
                else if (current_header_key == "SOURCE") append_opt(v.source_, val);
                else if (current_header_key == "ORGANISM") append_opt(v.organism_, val);
                else {
                    auto it = v.raw_fields_.find(current_header_key);
                    if (it != v.raw_fields_.end()) {
                        std::string t = trim(val);
                        if (!t.empty()) it->second += " " + t;
                    }
                }
            }

            continue;
        }

        if (section == Section::Features) {
            if (starts_with(line, "ORIGIN")) {
                if (in_qual_value) commit_multiline_qual();
                section = Section::Origin;
                continue;
            }

            // If we are inside a multiline quoted qualifier, append until closing quote.
            if (in_qual_value) {
                std::string cont = trim(slice_cols(line, 21));
                if (!cont.empty()) {
                    if (!q_value.empty()) q_value.push_back(' ');
                    q_value += cont;
                }

                if (q_quoted) {
                    // Close when we see a quote after the initial opening quote
                    std::size_t first = q_value.find('"');
                    std::size_t last  = q_value.rfind('"');
                    if (first != std::string::npos && last != std::string::npos && last > first) {
                        // Strip surrounding quotes
                        q_value = q_value.substr(first + 1, last - first - 1);
                        commit_multiline_qual();
                    }
                } else {
                    // Unquoted multiline is rare; commit now
                    commit_multiline_qual();
                }
                continue;
            }

            // Feature key at col 6..20 (index 5, len 15), location at col 22 (index 21)
            std::string fkey = trim(slice_cols(line, 5, 15));
            std::string rest = slice_cols(line, 21);

            if (!fkey.empty()) {
                Feature feat;
                feat.key = fkey;
                feat.location = trim(rest);
                v.features_.push_back(std::move(feat));
                current_feature = &v.features_.back();
                continue;
            }

            // Qualifiers start at col 22 and begin with '/'
            std::string qline = trim(rest);
            if (current_feature && !qline.empty() && qline[0] == '/') {
                qline.erase(0, 1); // remove leading '/'

                std::size_t eq = qline.find('=');
                if (eq == std::string::npos) {
                    // flag qualifier like /pseudo
                    current_feature->add_qualifier(qline, "");
                    continue;
                }

                q_key = qline.substr(0, eq);
                std::string rhs = trim(qline.substr(eq + 1));

                if (!rhs.empty() && rhs.front() == '"') {
                    q_quoted = true;
                    q_value = rhs;

                    std::size_t first = q_value.find('"');
                    std::size_t last  = q_value.rfind('"');
                    if (first != std::string::npos && last != std::string::npos && last > first) {
                        q_value = q_value.substr(first + 1, last - first - 1);
                        current_feature->add_qualifier(q_key, q_value);
                        q_key.clear();
                        q_value.clear();
                        q_quoted = false;
                    } else {
                        in_qual_value = true;
                    }
                } else {
                    q_quoted = false;
                    q_value = rhs;
                    current_feature->add_qualifier(q_key, q_value);
                    q_key.clear();
                    q_value.clear();
                }

                continue;
            }

            continue;
        }

        if (section == Section::Origin) {
            // Collect sequence letters only (normalize_sequence strips digits/whitespace anyway)
            v.sequence_ += normalize_sequence(line);
            continue;
        }
    }

    // Pick a name preference order for primer designer UX
    // (You can change this policy later.)
    if (v.definition_.has_value()) v.name_ = *v.definition_;
    else if (v.locus_.has_value()) v.name_ = *v.locus_;
    else if (v.accession_.has_value()) v.name_ = *v.accession_;
    else v.name_ = "GenBank_record";

    // Final cleanup (ensure normalized)
    v.sequence_ = normalize_sequence(v.sequence_);

    return v;
}

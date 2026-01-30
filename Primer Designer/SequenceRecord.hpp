#ifndef SEQUENCERECORD_HPP
#define SEQUENCERECORD_HPP

#include <cstddef>
#include <istream>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>


class SequenceRecord {
public:
    enum class Topology { Unknown, Linear, Circular };

    struct Feature {
        std::string key;        // e.g. "CDS", "gene", "primer_bind", "misc_feature"
        std::string location;   // raw GenBank location string, e.g. "935" or "complement(123..456)"
        std::unordered_map<std::string, std::vector<std::string>> qualifiers; // /note="...", /gene="..."

        // Helper: get first qualifier value if present
        std::optional<std::string> qualifier_first(const std::string& qkey) const;
        // Helper: add qualifier (supports repeated qualifiers)
        void add_qualifier(const std::string& qkey, const std::string& value);
    };

    // ---- Constructors ----
    SequenceRecord() = default;
    explicit SequenceRecord(std::string name, std::string sequence);

    // ---- Loading ----
    // Auto-detects based on content: FASTA if first non-ws char is '>', else tries GenBank.
    static SequenceRecord load_from_file(const std::string& path);

    static SequenceRecord load_fasta(std::istream& in);
    static SequenceRecord load_genbank(std::istream& in);

    // Convenience overloads by path
    static SequenceRecord load_fasta_file(const std::string& path);
    static SequenceRecord load_genbank_file(const std::string& path);

    // ---- Core accessors ----
    const std::string& name() const noexcept { return name_; }
    const std::string& sequence() const noexcept { return sequence_; }
    std::size_t length() const noexcept { return sequence_.size(); }
    Topology topology() const noexcept { return topology_; }

    void set_name(std::string n) { name_ = std::move(n); }
    void set_sequence(std::string seq); // normalizes (uppercases, strips spaces, converts U->T)
    void set_topology(Topology t) noexcept { topology_ = t; }

    // ---- Optional metadata (GenBank-ish) ----
    const std::optional<std::string>& locus() const noexcept { return locus_; }
    const std::optional<std::string>& definition() const noexcept { return definition_; }
    const std::optional<std::string>& accession() const noexcept { return accession_; }
    const std::optional<std::string>& version() const noexcept { return version_; }
    const std::optional<std::string>& source() const noexcept { return source_; }
    const std::optional<std::string>& organism() const noexcept { return organism_; }

    void set_locus(std::string v) { locus_ = std::move(v); }
    void set_definition(std::string v) { definition_ = std::move(v); }
    void set_accession(std::string v) { accession_ = std::move(v); }
    void set_version(std::string v) { version_ = std::move(v); }
    void set_source(std::string v) { source_ = std::move(v); }
    void set_organism(std::string v) { organism_ = std::move(v); }

    // Store any header fields not explicitly modelled (forward compatibility)
    const std::unordered_map<std::string, std::string>& raw_fields() const noexcept { return raw_fields_; }
    void set_raw_field(std::string key, std::string value);

    // ---- Features ----
    const std::vector<Feature>& features() const noexcept { return features_; }
    std::vector<Feature>& features_mut() noexcept { return features_; }

    void add_feature(Feature f);
    std::vector<const Feature*> find_features_by_key(const std::string& key) const;

    // ---- Primer-designer helpers ----
    // Returns true if sequence contains only A/C/G/T/N (after normalization); false if empty or has other chars.
    bool is_valid_dna() const;

private:
    // Core
    std::string name_;
    std::string sequence_;
    Topology topology_ = Topology::Unknown;

    // Optional metadata
    std::optional<std::string> locus_;
    std::optional<std::string> definition_;
    std::optional<std::string> accession_;
    std::optional<std::string> version_;
    std::optional<std::string> source_;
    std::optional<std::string> organism_;

    // Unknown/extra header fields
    std::unordered_map<std::string, std::string> raw_fields_;

    // Annotations
    std::vector<Feature> features_;

    // Internal helpers
    static std::string normalize_sequence(const std::string& s);
    static std::string trim(const std::string& s);
};

#endif // SEQUENCERECORD_HPP

# Primer Designer
Open-source C++ command-line tool for designing PCR primers for gene cloning, site-directed mutagenesis, and vector insertion, with automated primer quality checks and plasmid map graphic generation.  

## Description

primer-designer is a C++ program for designing PCR primers for common molecular biology workflows, including gene cloning, site-directed mutagenesis, and insertion of genes into plasmid vectors. The program operates entirely locally and parses standard FASTA sequence files along with JSON-formatted restriction enzyme libraries. Primer candidates are generated and filtered using standard design constraints, including length, melting temperature (calculated using both the Wallace rule and the SantaLucia (2004) nearest-neighbor model), GC content, GC clamps, and assessments of secondary structure such as hairpins and primer dimerization. For gene insertion workflows, the program supports restriction-enzyme-based strategies, allows user-specified tags, protease sites, linkers, and spacers, and produces linear restriction maps in text format and circular plasmid map graphics in SVG format. All functionality runs locally and completes in under half a second on modern hardware (e.g. Apple Silicon).

## Supported Workflows

### A) Generate flanking primers for gene cloning
**Input:**
- Gene sequence (as .fasta file in data/ directory)
  
**Output:**
- Text file containing a list of viable forward and reverse primers for cloning (cloning_primers<N>.txt)
---
### B) Generate mutagenic primers for site-directed mutagenesis by PCR
**Input:**
- Gene sequence (as .fasta file in data/ directory)
- Mutation specified in HGVS notation (e.g. c.120_121insT, c.86_91del, c.345C>T, c.100_110dup) (entered interactively via command line)

**Output:**
- Mutagenic primer pairs suitable for PCR-based mutagenesis (mutagenic_primers<N>.txt)
- FASTA file of the mutated gene sequence (mutated_gene<N>.fasta)
---
### C) Generate primers for insertion of a gene into a vector
**Input:**
- Gene sequence (as .fasta file in data/ directory)
- Vector sequence (as .fasta file in data/ directory)
- Restriction enzyme JSON file (already provided in data/RE_sites1.json)
  - If you would like to utilize a different list of restriction enzymes, use my repo restriction_enzyme_json to generate a new json file.
- User-specified design options (entered interactively via command line)
  - Choice of two restriction enzymes
  - Start and stop codons
  - N- and/or C-terminal tags (Built-in support: His6, FLAG, HA, Myc. Custom tags supported)
  - Protease cleavage sites (Built-in support: TEV, PreScission, Factor Xa, Enterokinase. Custom proteases supported)
  - Flexible linkers
  - Spacers (if necessary to maintain reading frame)

**Output:**
- Viable restriction enzyme combinations (viable_restriction_enzymes<N>.txt)
- Linear restriction map (linear_map<N>.txt)
- Circular plasmid map (circular_map<N>.svg)
- Primer pairs for gene insertion (gene_insertion_primers<N>.txt)
---
### D) Check for dimerization between two primers
 _* For use after generating and selecting viable primers using functions A-C._
 
**Input:**
- 2 primer sequences (entered interactively via command line)

**Output:**
- Answer to whether 2 primers heterodimerize (command line output)
---
### E) Generate the reverse complement of a sequence
 _* For general molecular biology usage._
 
**Input:**
- Nucleotide sequence (entered interactively via command line)

**Output:**
- Reverse complement of input sequence (command line output)


## Getting Started

### Requirements:
- A C++ compiler with support for modern C++ (e.g. Clang or GCC)
- macOS, Linux, or Windows (via a compatible toolchain)
- Standard C++ development environment
  (The program was developed in Xcode, but no platform-specific features are required.)

### Obtaining the source code:
Clone the repository from GitHub:
```
git clone https://github.com/chris-gud/primer-designer.git
cd primer-designer
```

### Building:
The program is written in standard C++ and can be built using Xcode or a command-line compiler.
Using Xcode (macOS):
- Open the project in Xcode.
- Select a build target.
- Build the project using the standard build command.
Using a command-line compiler (example):
```
g++ -std=c++17 -O2 -o primer-designer src/*.cpp
```
Compiler flags and source paths may be adjusted as needed for a given system.

### Executing program
The program is executed from the command line and operates interactively. Upon launch, the user is prompted to select a workflow and provide the required inputs, such as sequence file names, desired mutations, restriction enzymes, and optional design features. All .fasta and .json files must be placed in the data/ directory prior to use. 

Example input files are provided in the data/ directory:
```
data/
├── gene_example.fasta
├── vector_example.fasta
└── RE_sites1.json
```

All output files are written to an output/ directory, which is created automatically if it does not already exist.

Example output files are provided in the examples/ directory:
```
examples/
├── cloning/
│   └── cloning_primers_example.txt
├── mutagenesis/
│   ├── mutagenic_primers_example.txt
│   └── mutated_gene_example.fasta
└── gene_insertion/
    ├── viable_restriction_enzymes_example.txt
    ├── linear_map_example.txt
    ├── circular_map_example.svg
    └── gene_insertion_primers_example.txt
```

## Help



## Author
Chris Gudmundsen
2026

## Version History

* 1.0
    * Initial Release

## License

This project is licensed under the MIT License - see the LICENSE.md file for details

## Acknowledgments

* [List of Restriction Enzymes](https://www.neb.com/en-us/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities)
* [Example Gene (LacI)](https://www.ncbi.nlm.nih.gov/nuccore/671183183)
* [Example Vector (pAAV-AC-Myc-DDK AAV Expression Vector)](https://www.origene.com/catalog/vectors/aav-gene-expression-vectors/ps100089-paav-ac-myc-ddk-aav-gene-expression-vector)
* [SantaLucia 2004 Method](https://pubmed.ncbi.nlm.nih.gov/15139820/)
* [Mutation Notation](https://atlasgeneticsoncology.org/teaching/30067/nomenclature-for-the-description-of-mutations-and-other-sequence-variations#section-2)

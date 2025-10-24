# Nano Bisulfite Haplotag

A tool for haplotype tagging of bisulfite-converted nanopore sequencing data using SNP information.

## Features

- Support for bisulfite-converted nanopore sequencing reads
- Haplotype assignment based on heterozygous SNPs (hetSNPs)
- Configurable parameters for hetSNPs requirements
- Detailed statistics and logging
- Multi-processing support

## Installation

```bash
git clone https://github.com/sherryxuePKU/nano-bisulfite-haplotag.git
cd nano-bisulfite-haplotag
pip install -e .
```

## Input File Formats

> BAM File Format  

The input BAM file must meet the following requirements:

- Coordinate sorted: Use samtools sort to sort by coordinates
- Indexed: BAM index file (.bai) must be present in the same directory
- Bisulfite-converted reads: Reads should be from bisulfite sequencing data
- Proper alignment: Reads must be aligned to reference genome  

> SNP File Format

The SNP file requires following format with tab-separated values: 

**Example SNP file**

```
chr1    1000    A    G    rs123456
chr1    2500    C    T    rs234567
chr1    5000    G    A    rs345678
chr2    1200    T    C    rs456789
chr2    3400    A    T    rs567890
```

**Column descriptions**
- Column 1 (CHROMOSOME): Chromosome name (e.g., chr1, chr2, 1, 2)
- Column 2 (POSITION): 1-based genomic position
- Column 3 (REF_ALLELE): Reference allele (A, T, G, C)
- Column 4 (ALT_ALLELE): Alternative allele (A, T, G, C)
- Column 5+ (OPTIONAL): Additional annotation columns (SNP ID, etc.)

## Usage

### Common line

```
nano-bisulfite-haplotag -b input.bam -s snps.txt -o output_dir
```

### Python API

```
from nano_bisulfite_haplotag import HaploTagger

tagger = HaploTagger(
    bam_file="input.bam",
    snp_file="snps.txt",
    output_dir="output",
    min_snp=2,
    hap_fold_change=2.0,
    num_processors=4
)

stats = tagger.analyze_snps()
print(f"Tagging rate: {stats['tagging_rate']:.2%}")
```

### Parameters

* `-b, --input_bam`: Input BAM file (coordinate sorted with index)
* `-s, --snp_file`: SNP file in SNPsplit format
* `-o, --output_dir`: Output directory (default: BAM file location)
* `-m, --min_snp`: Minimum SNPs per read (default: 2)
* `-f, --hap_fc`: Minimum fold change for haplotype assignment (default: 2.0)
* `-n, --num_processors`: Number of processors (default: 4)

## Output Files

* `*_haplotag_results.tsv`: Complete results with haplotype assignments
* `*_HAP1_reads.tsv`: Reads assigned to haplotype 1
* `*_HAP2_reads.tsv`: Reads assigned to haplotype 2
* `*_haplotag_stats.txt`: Analysis statistics
* `haplotag_*.log`: Detailed log file

**Example Results file**

```
read_id         chromosome    start     end       haplotype    hap1_snps    hap2_snps    total_snps    confidence
read_001        chr1          1000      2000      HAP1         5            0            5             1.0
read_002        chr1          1500      2500      HAP2         0            3            3             1.0
read_003        chr1          2000      3000      AMBIGUOUS    2            2            4             0.5
```

## Requirements

* Python ≥3.7
* pysam ≥0.19.0
* pandas ≥1.3.0
* numpy ≥1.20.0

## Citation

If you use this tool in your research, please cite:

```
Lin J, Xue X, Wang Y, et al. scNanoCOOL-seq: a long-read single-cell sequencing method for multi-omics profiling within individual cells. Cell Res. 2023;33(11):879-882. doi:10.1038/s41422-023-00873-5
```

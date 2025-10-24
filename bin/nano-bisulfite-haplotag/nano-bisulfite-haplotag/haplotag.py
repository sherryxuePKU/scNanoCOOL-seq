"""
Core haplotype tagging functionality
"""

import os
import logging
import pysam
from typing import Dict, List, Tuple, Optional
from multiprocessing import Pool
import pandas as pd


class HaploTagger:
    """
    Main class for performing haplotype tagging on bisulfite sequencing data.
    """
    
    def __init__(self, 
                 bam_file: str, 
                 snp_file: str, 
                 output_dir: str,
                 min_snp: int = 2,
                 hap_fold_change: float = 2.0,
                 num_processors: int = 4):
        """
        Initialize HaploTagger.
        
        Args:
            bam_file: Path to input BAM file (coordinate sorted with index)
            snp_file: Path to SNP file (SNPsplit format)
            output_dir: Output directory for results
            min_snp: Minimum SNPs covered by single reads
            hap_fold_change: Minimum fold change of haplotypes
            num_processors: Number of processors for analysis
        """
        self.bam_file = bam_file
        self.snp_file = snp_file
        self.output_dir = output_dir
        self.min_snp = min_snp
        self.hap_fold_change = hap_fold_change
        self.num_processors = num_processors
        
        # Validate inputs
        self._validate_inputs()
        
        # Initialize logging
        self.logger = logging.getLogger(__name__)
        
        # Load SNPs
        self.snps = self._load_snps()
        
    def _validate_inputs(self) -> None:
        """Validate input files and parameters."""
        if not os.path.exists(self.bam_file):
            raise FileNotFoundError(f"BAM file not found: {self.bam_file}")
        
        if not os.path.exists(self.snp_file):
            raise FileNotFoundError(f"SNP file not found: {self.snp_file}")
        
        if not os.path.exists(f"{self.bam_file}.bai"):
            raise FileNotFoundError(f"BAM index not found: {self.bam_file}.bai")
        
        if self.min_snp < 1:
            raise ValueError("min_snp must be >= 1")
        
        if self.hap_fold_change < 1.0:
            raise ValueError("hap_fold_change must be >= 1.0")
        
        if self.num_processors < 1:
            raise ValueError("num_processors must be >= 1")
    
    def _load_snps(self) -> Dict[str, List[Tuple[int, str, str]]]:
        """
        Load SNPs from file.
        
        Returns:
            Dictionary with chromosome as key and list of (pos, ref, alt) as values
        """
        self.logger.info(f"Loading SNPs from {self.snp_file}")
        
        snps = {}
        with open(self.snp_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                
                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[2]
                alt = parts[3]
                
                if chrom not in snps:
                    snps[chrom] = []
                
                snps[chrom].append((pos, ref, alt))
        
        # Sort SNPs by position
        for chrom in snps:
            snps[chrom].sort(key=lambda x: x[0])
        
        total_snps = sum(len(snps[chrom]) for chrom in snps)
        self.logger.info(f"Loaded {total_snps} SNPs across {len(snps)} chromosomes")
        
        return snps
    
    def _process_read(self, read: pysam.AlignedSegment) -> Optional[Dict]:
        """
        Process a single read to determine haplotype assignment.
        
        Args:
            read: pysam AlignedSegment object
            
        Returns:
            Dictionary with read information and haplotype assignment, or None
        """
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            return None
        
        chrom = read.reference_name
        if chrom not in self.snps:
            return None
        
        read_start = read.reference_start
        read_end = read.reference_end
        
        # Find overlapping SNPs
        overlapping_snps = []
        for pos, ref, alt in self.snps[chrom]:
            if read_start <= pos <= read_end:
                overlapping_snps.append((pos, ref, alt))
        
        if len(overlapping_snps) < self.min_snp:
            return None
        
        # Count haplotype-supporting SNPs
        hap1_count = 0  # Reference allele
        hap2_count = 0  # Alternative allele
        
        for pos, ref, alt in overlapping_snps:
            # Get base at SNP position
            read_pos = read.query_position_or_next(pos - read.reference_start)
            if read_pos is None:
                continue
            
            read_base = read.query_sequence[read_pos].upper()
            
            if read_base == ref.upper():
                hap1_count += 1
            elif read_base == alt.upper():
                hap2_count += 1
        
        total_snps = hap1_count + hap2_count
        if total_snps < self.min_snp:
            return None
        
        # Determine haplotype assignment
        if hap1_count == 0 and hap2_count > 0:
            haplotype = "HAP2"
            confidence = 1.0
        elif hap2_count == 0 and hap1_count > 0:
            haplotype = "HAP1"
            confidence = 1.0
        elif hap1_count > 0 and hap2_count > 0:
            fold_change = max(hap1_count, hap2_count) / min(hap1_count, hap2_count)
            if fold_change >= self.hap_fold_change:
                haplotype = "HAP1" if hap1_count > hap2_count else "HAP2"
                confidence = fold_change / (fold_change + 1)
            else:
                haplotype = "AMBIGUOUS"
                confidence = 0.5
        else:
            return None
        
        return {
            'read_id': read.query_name,
            'chromosome': chrom,
            'start': read_start,
            'end': read_end,
            'haplotype': haplotype,
            'hap1_snps': hap1_count,
            'hap2_snps': hap2_count,
            'total_snps': total_snps,
            'confidence': confidence
        }
    
    def analyze_snps(self) -> Dict[str, int]:
        """
        Perform SNP analysis and haplotype tagging.
        
        Returns:
            Dictionary with analysis statistics
        """
        self.logger.info("Starting haplotype analysis")
        
        results = []
        processed_reads = 0
        tagged_reads = 0
        
        # Open BAM file
        with pysam.AlignmentFile(self.bam_file, 'rb') as bam:
            for read in bam.fetch():
                processed_reads += 1
                
                if processed_reads % 10000 == 0:
                    self.logger.info(f"Processed {processed_reads} reads, tagged {tagged_reads}")
                
                result = self._process_read(read)
                if result is not None:
                    results.append(result)
                    tagged_reads += 1
        
        self.logger.info(f"Analysis complete. Processed {processed_reads} reads, tagged {tagged_reads}")
        
        # Save results
        self._save_results(results)
        
        # Generate statistics
        stats = self._generate_statistics(results, processed_reads)
        self._save_statistics(stats)
        
        return stats
    
    def _save_results(self, results: List[Dict]) -> None:
        """Save haplotype tagging results."""
        if not results:
            self.logger.warning("No reads were successfully tagged")
            return
        
        # Convert to DataFrame
        df = pd.DataFrame(results)
        
        # Save main results
        output_file = os.path.join(self.output_dir, f"{os.path.basename(self.bam_file)}_haplotag_results.tsv")
        df.to_csv(output_file, sep='\t', index=False)
        self.logger.info(f"Results saved to {output_file}")
        
        # Save haplotype-specific files
        for hap in ['HAP1', 'HAP2']:
            hap_df = df[df['haplotype'] == hap]
            if not hap_df.empty:
                hap_file = os.path.join(self.output_dir, f"{os.path.basename(self.bam_file)}_{hap}_reads.tsv")
                hap_df.to_csv(hap_file, sep='\t', index=False)
                self.logger.info(f"{hap} reads saved to {hap_file}")
    
    def _generate_statistics(self, results: List[Dict], total_reads: int) -> Dict[str, int]:
        """Generate analysis statistics."""
        if not results:
            return {
                'total_reads': total_reads,
                'tagged_reads': 0,
                'hap1_reads': 0,
                'hap2_reads': 0,
                'ambiguous_reads': 0,
                'tagging_rate': 0.0
            }
        
        df = pd.DataFrame(results)
        
        stats = {
            'total_reads': total_reads,
            'tagged_reads': len(results),
            'hap1_reads': len(df[df['haplotype'] == 'HAP1']),
            'hap2_reads': len(df[df['haplotype'] == 'HAP2']),
            'ambiguous_reads': len(df[df['haplotype'] == 'AMBIGUOUS']),
            'tagging_rate': len(results) / total_reads if total_reads > 0 else 0.0
        }
        
        return stats
    
    def _save_statistics(self, stats: Dict[str, int]) -> None:
        """Save analysis statistics."""
        stats_file = os.path.join(self.output_dir, f"{os.path.basename(self.bam_file)}_haplotag_stats.txt")
        
        with open(stats_file, 'w') as f:
            f.write("Haplotype Tagging Statistics\n")
            f.write("============================\n\n")
            
            for key, value in stats.items():
                if key == 'tagging_rate':
                    f.write(f"{key}: {value:.4f}\n")
                else:
                    f.write(f"{key}: {value}\n")
        
        self.logger.info(f"Statistics saved to {stats_file}")
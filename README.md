# Split10x

Given a bam file from 10x genomics with barcode inforamtion 
and a phased-vcf (again by 10x), generate two bams that contain 
the reads from each haplotype.

Usage:
split10x reads.bam vars.vcf hap_base max_molecule_length


max_molecule_length is the maximum expected size of a barcoded fragment of DNA.

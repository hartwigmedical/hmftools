
# Somatic Alterations in Genome (SAGE)

SAGE looks directly at the provided reference and tumor bams looking for evidence of:
1. Inframe indels within the specified coding regions; and
2. Specific hotspots at known locations.

All sites examined are written to an output VCF with some soft filtering applied to the results.  

## Sites
The supplied set of known hotspots is combined with any inframe indels found in the tumor to form the complete set of sites for which to collect evidence. 

To find inframe indel sites, SAGE examines the CIGAR field of all tumor alignments within the supplied coding regions looking for inserts or deletes that are divisible by 3, e.g., 100M3D50M and 60M9I81M.

## Accumulating Evidence
In an alignment, evidence of the alt is accumulated only if:
1. Quality of the read is sufficient; and
2. The alt matches exactly. 

The quality of an insert, SNV or MNV is calculated as the average quality of all alt bases. The quality of a delete is taken from the base after the delete if available, otherwise from the base prior. The minimum acceptable quality is set with min_base_quality [13].

SNVs and MNVs cannot be part of a larger MNV, eg, C > T does not match CA > TT. There must be at least one base either side of the variant that match the reference exactly. 

Only simple indels are considered (eg A > AC). Evidence of complex indels (eg, A > CC) will not accumulate. 

## Filtering Evidence
When writing to a file, a number of soft filters are applied. 

Any variant will be filtered as **LOW_CONFIDENCE** under the following conditions:
1. Insufficient tumor allelic depth, i.e., tumor allelic depth < min_tumor_reads [2]; or
2. Evidence of support for the ALT in the germline sample, ie. filter if germline ALT read count >1 OR if germline ALT read count = 1 AND germline heterozygous binomial likelihood > max_het_binomial_likelihood [0.01]

Any snv or mnv will be filtered as **LOW_CONFIDENCE** if:
1. Insufficient quality, i.e., quality < min_snv_quality [100]; or
2. Insufficient VAF, i.e., vaf < min_snv_vaf [0.005]

Any indel will be filtered as **LOW_CONFIDENCE** if:
1. Insufficient quality, i.e., quality < min_indel_quality [150]; or
2. Insufficient VAF, i.e., vaf < min_indel_vaf [0.02]

Any unfiltered indel (regardless of being a hotspot or not) will be filtered as **GERMLINE_INDEL** if there exists an indels at that site in the germline. 

## Output VCF
Only one variant per site is written to file as ranked by filter then quality.


## Post Processing
The HMF pipeline employs some post processing steps to further filter the output. We have generated a PON file from 1077 sage results that captures any variants that have more than one read in the germline of at least 2 samples. These variants are filtered from our output.

## Usage

Argument | Description 
---|---
coding_regions | Coding regions bed file to search for inframe indels
known_hotspots | Tab separated file of known hotspot locations
max_het_binomial_likelihood | Up to 1 read of support for the ALT in the reference is allowed if the binomial likelihood of observing the read (assuming a 50% VAF) is less than this parameter [0.01]
min_base_quality | Minimum quality for a base to be considered [13]
min_indel_quality | Low confidence filtering minimum indel quality [150]
min_indel_vaf | Low confidence filtering minimum indel VAF [0.02]
min_mapping_quality |  Minimum mapping quality for an alignment to be used [1]
min_snv_quality | Low confidence filtering minimum SNV/MNV quality [100]
min_snv_vaf | Low confidence filtering minimum SNV/MNV VAF [0.005]
min_tumor_reads | Low confidence filtering minimum tumor reads [2]
out | Output VCF file to write. Gz  supported
ref_genome | Path to the ref genome fasta file
reference | Name of reference sample
reference_bam | Path to reference bam file
tumor | Name of tumor sample
tumor_bam | Path to tumor bam file
 

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

SNVs and MNVs cannot be part of a larger MNV, eg, C > T does not match CA > TT. There must be at least two bases either side of the variant that match the reference exactly. 

Only simple indels are considered (eg A > AC). Evidence of complex indels (eg, A > CC) will not accumulate. 

## Collating Evidence
Although multiple variants may start at the same position, only one will be written to file. Variants are ranked first by size (favouring length) then by quality. 

## Filtering Evidence
When writing to a file, a number of soft filters are applied. 

Any variant will be filtered as **LOW_CONFIDENCE** under the following conditions:
1. Insufficient tumor allelic depth, i.e., tumor allelic depth < min_tumor_reads [2]; or
2. Excessive germline allelic depth, i.e., germline allelic depth > 2 OR germline allelic depth = 1 AND germline heterozygous binomial likelihood > max_het_binomial_likelihood [0.01]

Any hotspot will be filtered as **LOW_CONFIDENCE** if:
1. Insufficient quality, i.e., quality < min_hotspot_quality [100]; or
2. Insufficient VAF, i.e., vaf < min_hotspot_vaf [0.005]

Any inframe indel (that is not also a hotspot) will be filtered as **LOW_CONFIDENCE** if:
1. Insufficient quality, i.e., quality < min_inframe_quality [150]; or
2. Insufficient VAF, i.e., vaf < min_inframe_vaf [0.02]

Any unfiltered indel (regardless of being a hotspot or not) will be filtered as **GERMLINE_INDEL** if there exists an indels at that site in the germline. 

## Post Processing
The HMF pipeline employs some post processing steps to further filter the output. We have generated a PON file from 1077 sage results that captures any variants that have more than one read in the germline of at least 2 samples. These variants are filtered from our output.

 

# Somatic Alterations in Genome (SAGE)
SAGE is a somatic SNV, MNV and small INDEL caller optimised to call narrow regions of the genome with high prior chance of a variant with very high sensitivity. SAGE looks directly at the provided reference and tumor bams looking for evidence of:
1. Inframe indels within the specified coding regions; and
2. Specific hotspots at known locations.

All sites examined are written to VCF with some soft filtering applied to the results.  

## Sites
The supplied set of known hotspots is combined with any inframe indels found within specified regions of the tumor to form the complete set of sites for which to collect evidence. 

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

Any unfiltered indel will be filtered as **GERMLINE_INDEL** if there exists an indel at that site in the germline. 

## Output VCF
Only one variant per site is written to file as ranked by filter then quality.

A snippet of output follows:

```
##fileformat=VCFv4.2
##FILTER=<ID=GERMLINE_INDEL,Description="Set if indel has any germline indels at that site.">
##FILTER=<ID=LOW_CONFIDENCE,Description="Set if excessive germline reads or insufficient quality or tumor reads">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allelic Frequency">
##INFO=<ID=GHBL,Number=1,Type=Float,Description="Germline Heterozygous Binomial Likelihood. Only applies when germline alt support = 1">
##INFO=<ID=HOTSPOT,Number=1,Type=String,Description="Hotspot Type: known, inframe">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	COLO829R	COLO829T
1	9780851	.	GAG	AAA	0	LOW_CONFIDENCE	AF=0.00;HOTSPOT=known	GT:AD:DP	0/1:31,0:32	0/1:56,0:58
5	1295228	.	GG	AA	1113	PASS	AF=0.917;HOTSPOT=known	GT:AD:DP	0/1:17,0:20	0/1:0,33:36
7	128829039	.	GGCT	G	203	PASS	AF=0.055;HOTSPOT=inframe	GT:AD:DP	0/1:36,0:36	0/1:103,6:109
7	140453136	.	A	T	4047	PASS	AF=0.642;HOTSPOT=known	GT:AD:DP	0/1:35,0:35	0/1:56,104:162
11	1018167	.	ATCG	A	253	GERMLINE_INDEL	AF=0.051;GHBL=0.00;HOTSPOT=inframe	GT:AD:DP	0/1:71,1:72	0/1:148,8:156
```

## Arguments

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

##Example Usage
```
java -Xmx8G -Xms4G \
    -cp sage.jar com.hartwig.hmftools.sage.SageHotspotApplication \
    -tumor COLO829T -tumor_bam /path/to/COLO829T.bam \
    -reference COLO829R -reference_bam /path/to/COLO829R.bam \
    -known_hotspots /path/to/KnownHotspots.tsv \
    -coding_regions /path/to/CodingRegions.bed \
    -ref_genome /path/to/refGenome.fasta \
    -out COLO829T.vcf.gz
```

## Inputs

Hotspots file should be header-less and tab separated with columns: Chromosome, Position, Ref, Alt, ie:
``` 
1	8073509	C	A
1	9777666	C	A
1	9777666	C	G
1	9780851	G	A
``` 

## Post Processing
The HMF pipeline employs some post processing steps to further filter the output. We have generated a PON file from 1077 sage results that captures any variants that have more than one read in the germline of at least 2 samples. These variants are filtered from our output.

 
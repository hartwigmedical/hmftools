# WISP

WISP estimates TF in a given sample (typically ctDNA), guided by the biomarkers identified prior anlaysis of the same patient (typically a primary tissue biopsy).   The algorithm supports data generated from both WGS or tarrgeted analyses.    

## Commands

```
<TO DO>
```

## Arguments

Argument | Required | Description
---|---|---
<TO DO>| |

## Algorithm

WISP can infer tumor fraction of a sample via 3 different meths
- SNV (& small Indel)
- CNA (WGS only)
- LOH (WGS only)

###  SNV tumor fraction estimate
Using the cfDNA annotated tissue biopsy purple VCF, apply the following filters 
- High confidence variant in primary:  Mappability>=0.5 & RepeatCount<4 
- High confidence observation in ctDNA: QualPerAD>minQualPerAD|AD==0 
- Not subclonal:  AND(SubclonalPerc<0.5,VCN<0.7) 
- GC content:  GC[probe] > 0.4 (if in PANEL mode) 

Considering only the remaining unfiltered sites calculate the TFctDNA using the following formula 
```
adjVAFctDNA = ADctDNA/DPctDNA -ε 
Where ε = expected noise in cfDNA  
TFctDNA = 2* adjVAFctDNA / (WA_VCN + 2* adjVAFctDNA - WA_CN * adjVAFctDNA) 
```

The expected noise is set as the depth weighted average empirical error rate from SAGE BQR of all fragments based on the trinucleotide context for variants with high raw base qual (>=30).    

A check is made at this stage for variants which are outliers and may unduly contribute to the fit.   Variants with a VAF of more than 8x the average cfDNA VAF for the sample, with allele support of at least 5 and whose allele support makes up more than 10% of reads in the sample are filtered as outliers.  These are excluded as likely germline variants which have leaked into the somatic variants or clonal haematopoiesis. The ctDNA_TF is then recalculated excluding the filtered observation. 
 
The estimated ctDNA_TF represents the purity if all SNV considered from the tissue biopsy are still present in the ctDNA. However, the % of SNV that are present is frequently observed to be much lower than this, since the tissue biopsy may be taken from a locally dominant subclone (eg. a single metastasis) which has a number of private mutations.  The dominant clone present in the cfDNA may have evolved from an earlier or different subclone in the primary. Hence the raw ctDNA_TF at this stage represents a lower bound for the purity estimate.  

One of 2 methods are used to estimate the proportion of variants which are clonal depending on the depth and allelic depth 

#### VAF_PEAK clonality method
If the weightedAverageDepth (wAD^) of variants in the sample is greater than 20 AND at least 10 points have more than 1 fragment support AND there are EITHER [at least 3 times as many sites with more than 1 support than 1 fragment support (ie. [N>1]/N[1] >3) or rawTF>1% and [N>1]/N[1] >3], then a kernel density estimator is used to find a TF peak which may be higher than the TF. Only variants with DP > max(20, 10% of wAD) are used. For each individual variant, an implied TF is calculated.  Variants with very high implied TF are capped at 200%.  The bandwidth is set to:   
```
bw = min[max(0.2*rawTF,5/wAD, 0.125 * nth highest TF), 3*rawTF,0.1] 
Where n = max (4,8th percentile of variants) 
```
The final TF is set to the highest TF peak with weight > max(3,3%) variants  
#### LOW_COUNT  clonality method
If the wAD < 50,  and there are at least max(4,0.5% of variants) with 2 or more fragment support, then the ratio of site count with 1 fragment support (N[1]) compared to the site count with >1 fragment support (N[>1]) is used to estimate the support.  Only variants with VCN <= 2* median VCN and sampleDP <=2*wAD are used.   Specifically, the model assumes that a proportion of variants has been dropped from the sample (dropout%) and tests a matrix of combinations of {sampleVAF, dropout%} to see which best explains the ratio N[>1]/N[1].   The combination which best matches the observed ratio is selected.  

^Note that, wAD is weighted for CN and depth in the tissue sample and is defined as  
```
wAD = Σ(i=1->n)[(DPi_cfDNA)2/CNn] / Σ(i=1->n)[DPi_cfDNA/CNi] 
```
The adjusted VAF is used to recalculate the implied tumor fraction estimate for the sample.  If the conditions for neither model are met, then the VAF is not adjusted.    

#### Other annotations
- **Tumor Presence P-Value** – Poisson test of whether the observed AD is greater than expected noise.  Since some trinucleotides can have much higher rates than others (particularly C>T), WISP tests for tumor presence by filtering sites based on the empirically measured error rates in SAGE BQR with raw base qualities higher than 35, for each context at 4 different phred score thresholds f {All sites, >38,>42,>45}.  For each threshold, the observed AD is compared the weighted average expected error rate for the included sites (note the maximum BQR phred score is limited to 50 per variant context to be able to handle sparse data).   The final p-value is set to the minimum threshold across all the tests.   For duplex sequencing a further test is made to see if the DUAL_STRAND AD > expected noise (set to a fixed 1e-6 for DUPLEX observations). 
- **LOD** – This is calculated as the observed TF that would return a 99% confidence of tumor presence (ie p-value = 0.01) given the allele fragments. Again, the lowest LOD across different phred thresholds is used. 
- **LowTFEstimate** – Estimated using the same formula with the following adjustments: 
   - Recalculate TF based on AD = poisson[sumAD,0.05]  
   - VAF_PEAK : double the bandwidth and multiply by lowRawTFEstimate/rawTF 
   - LOW_COUNT : use max(0,min[poisson(N[>1],0.05), N[>1]-2] in lieu of N[>1] algorithm and multiply by lowRawTFEstimate/rawTF 
- **HighTFEstimate**: Estimated using the same formula with the following adjustments: 
   - Recalculate TF based on AD = poisson[sumAD,0.95]  
   - VAF_PEAK: halve the bandwidth and multiply by highRawTFEstimate/rawTF 
   - LOW_COUNT : use max[poisson(N[>1],0.95), N[>1]+2] in lieu of N[>1] and multiply by highRawTFEstimate/rawTF 

### CNA tumor fraction estimate (WGS only) 
 
The tumor fraction is also estimated independently using a somatic CNA approach, by comparing the median tumorGCratio observed in cobalt in the cfDNA to the estimated copy number for each distinct copy number segment identified in the tumor biopsy.    Note that the CNA based estimate is only run for WGS cfDNA. 

To do this, first we filter for segments with integer copy numbers (ie within 0.2 of an integer value) between copy number 0 and 6.   Non integer copy numbers likely indicate sub-clonal copy numbers which may be volatile over the evolution of the tumor and hence are excluded from consideration.  High copy number regions may also be volatile in copy number (eg. EcDNA) and are also excluded).    

For each distinct copy number segment with these characteristics, the tumorGCratio is calculated for the cfDNA WGS sample (using all depth windows with mappability > 0.85 and with 0.2<GCContent<0.6), and the count of depth windows is used to weigh the observations. For each integer CN level in tissue biopsy (0 to 6), a weighted median of the tumorGCratio in the cfDNA is calculated across all the segments sharing that copy number level.   A regression is then performed to find the purity of the cfDNA that best fits the observed tumorGCRatio pattern. A confidence interval is determined by rerunning the same regression but excluding each chromosome one at a time.  The highest and lowest observations are set to the confidence interval for the purity estimate 

Note that if the copy number variation in the tissue biopsy may be minimal.   To deal with this we measure an aneuploidyScore = 1 - (max % of clonal depth windows  at single CN level).  If aneuploidyScore<0.001 the copy number purity is set to 0.    We found empirically using unrelated controls that the CN approach is subject to noise if the anneuploidy is low OR if a high proportion of the sample is subclonal (which may indicate either a poor fit or very unstable primary tumor). We therefore estimate an approximate limit of detection based on these 2 factors using the following model: 
```
LOD_CNA = 0.4% / sqrt(aneuploidy score) / clonalProportion^3 
```
### LOH estimate (WGS only) 
In WGS mode, a tumor fraction estimate is also made based on the LOH sites.  WISP considers each high confidence heterozygous BAF site in an LOH region identified in the primary tumor by PURPLE.   If the lost allele can be determined with high confidence (tumorBAF < 0.45 => ALT lost ; tumorBAF>0.55 => REF lost), then the site is used in the LOH estimate.   

An AF is first calculated per region using the high confidence LOH sites as  
```
AF = ΣAD[LostAllele]/ ΣDP  
```
If tumor is present, then the AF of the lost allele foreach region is expected to be below 50%.   Regions with average AF > 0.55 are filtered from subsequent calculations as the allelic balance is clearly not in the correct direction.  

The implied tumor fraction for each region is estimated as: 
```
TF = (1-2*AF) / [AF*(CN-2)+1] 
```
WISP then uses the median TF across all the regions as the final TF (weighted by number of sites per region).   A final check is made to see if there is a subset of regions which have TF significantly higher than the median indicating that some of the identified LOH regions may have been subclonal in the tumor.  For this WISP searches for regions with TF > max(median*1.2, median + 2%) AND impliedTF >> medianTF (p-value).   If at least 1% of the genome satisifies these condtions, then the TF is instead set to the median of those selected regions.    

A LOD is calculated based on the TF consistent with an AF that would return a p-value of <0.01 given the total fragments observed at LOH sites.

## Output Files
### Summary output
Field | Description 
---|---
PatientId | Id of Patient 
SampleId | Id of sample 
TumorPurity | Purity of primary tissue biopsy 
TumorPloidy | Ploidy of primary tissue biopsy
TotalVariants | Total number of variants evaluated for SNV TF estimate 
CalcVariants | # sites passing all filters for SNV TF estimate 
ChipVariants | # of sites excluded as clonal hematopeioisis variants 
SNVPurity | Final SNV cfDNA tumor fraction estimate (adjusted for clonality) 
RawSNVPurity | Raw SNV cfDNA tumor fraction estimate (not adjusted for clonality) 
SNVPValue | P-value likelihood of tumor presence in sample based on SNV support 
SNVPurityLow | Lower bound of SNV purity estimate 
SNVPurityHigh | Upper bound of SNV purity estimate 
ClonalMethod | Clonality method used = {VAF_PEAK,LOW_COUNT, NONE} 
Frag1Variants  | # sites passing all filters with 1 fragment support 
Frag2PlusVariants | # sites passing all filters with >1 fragment support 
ClonalPeakVariants | # sites found in clonal peak (VAF_PEAK clonality method) 
ClonalDropoutRate | % of sites estimated to be dropped (LOW_COUNT clonality method) 
DualSNVPValue | P-value likelihood of tumor presence in sample based on DUAL_STRAND SNV support (may be lower than global estimate due to lower error rates for DUAL_STRAND) 
SNVLod | Estimated minimum tumor fraction limit of detection based on number of fragments overlapping included SNV sites 
TotalFragments  | Total depth of coverage across all included SNV sites 
TotalDual  | Depth of DUAL_STRAND coverage across all included SNV sites 
AlleleFragments | Fragments supporting ALT across all included SNV sites 
AlleleDual  | DUAL_STRAND fragments supporting ALT across all included SNV sites 
WeightedAvgDepth | Average depth weigted by AD support 
WeightedAvgVCN | Weighted average variant copy number 
WeightedAvgCN  | Weighted average copy number 
PeakBandwidth  | Bandwidth used for VAF_PEAK clonality method 
PeakBandwidthLow | Lower bound bandwidth used for VAF_PEAK clonality method 
PeakBandwidthHigh | Upper bound bandwidth used for VAF_PEAK clonality method 
ErrorRate | Error rate used in SNV TF estimate (weighted average of included variants), based on the weighted average error rate of sites selected 
RawBQRErrorRate | Average empirical error rate for the sample. 
BqrThreshold | Minimum raw base quality used to filter .sites (either 0,38,42,45 is chosen based on whichever has the lowest LOD) 
BQRExtraInfo | Error rates and calculated LOD per BQR filter threshold. 
LOHPurity | Estimated purity from LOH 
LOHSiteCount | Total number of sites used in LOH algorithm 
LOHPercent | Percentage of primary sample BAF points with LOH 
LOHMeanCopyNumber | Mean copy number across all LOH BAF points 
LOHMedianAF | Weighted median AF across LOH points 
LOHMeanAF | MeanAF across LOH poiints 
LOHPValue  | Probability of tumor presence based on LOH algorithm 
LOHFragments | Total allele fragments used in LOH algorithm 
CopyNumberPurity  | Estimate of cfDNA tumor fraction from copy number (CNA) 
CopyNumberPurityLow  | Lower bound of CNA cfDNA tumor fraction estimate 
CopyNumberPurityHigh | Upper bound of CNA cfDNA tumor fraction estimate 
AnueploidyScore | Proportion of tissue biopsy genome with CN != most frequent CN 
ClonalPercent | Proportion of tissue biopsy genome with CN not with +/- 0.2 of an integer value 

### Somatic Variant 
This file contains counts and filtering rules for each variant used in the SNV tumor fraction estimate 

Field | Description 
---|---
PatientId | Id of Patient 
SampleId | Id of sample 
Chromosome | Chromosome of variant 
Position | Position of variant 
Ref | Ref nucleotide of variant 
Alt | Alt nucleotide of variant 
IsProbe | Variant was included in panel design 
Filter | Either PASS or filter reason 
Tier | SAGE variant calling tier {HOTSPOT,PANEL,HIGH_CONFIDENCE,LOW_CONFIDENCE} 
Type | {INDEL,MNV or SNV} 
TNC | Trinucleotide context of variant (used for error rate estimation) 
Mappability | Mappability of variant position 
SubclonalPerc | Probability that variant is subclonal 
RepeatCount | Repeat count of variant 
Gene | Gene context of variant 
CodingEffect | {MISSENSE,NONSENSE_OR_FRAMESHIFT,SPLICE or NONE} 
Hotspot | Is hotspot {TRUE,FALSE} 
Reported | Variant is reportable {TRUE,FALSE} 
VCN | Variant copy number (from primary tumor) 
CopyNumber | Copy number of variant position (from primary tumor) 
TumorDP | Depth at variant site in primary tumor 
TumorAD | Allelic depth support for variant in primary tumor 
SampleDP | Depth at variant site in cfDNA sample 
SampleAD | Allelic depth support for variant in cfDNA sample 
SampleDualDP | DUPLEX depth at variant site in cfDNA sample 
SampleDualAD | DUPLEX allelic depth support for variant in cfDNA sample 
QualPerAD | Average quality per allelic depth support in SAGE (used as filter) 
SeqGCRatio | GC ratio of 120 bases surrounding variant site 
BQRErrorRate | Empirically measure error rate for variant trinucleotide + alt context from SAGE. 

### SCNA output 
Field | Description 
---|---
PatientId | Id of Patient 
SampleId | Id of Sample 
Chromosome | Chromosome of copy number segment in primary tumor 
SegmentStart | Start position of segment 
SegmentEnd | End position of segment 
CopyNumber | Copy number of segment 
GcRatioCount | Depth windows in segment in primary tumor 
GcRatioMedian | Median GC ratio of segment in cfDNA sample 
GcRatioMean | Mean GC ratio of segment in cfDNA sample 




## Known Issues / Future improvements
Additional biomarkers
- SV
- C-circles
- CDR3

Additional features
- Fragmentomics: 5' & 3' trinculeotide frequencies AND positional biases of ref and alt fragments per site
- TINC: add explicit support
- Support for phased analysis 


 ## Version History and Download Links


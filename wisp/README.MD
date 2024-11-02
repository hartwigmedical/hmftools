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

#### VAF_PEAK
If the weightedAverageDepth (wAD^) of variants in the sample is greater than 20 AND at least 10 points have more than 1 fragment support AND there are EITHER [at least 3 times as many sites with more than 1 support than 1 fragment support (ie. [N>1]/N[1] >3) or rawTF>1% and [N>1]/N[1] >3], then a kernel density estimator is used to find a TF peak which may be higher than the TF. Only variants with DP > max(20, 10% of wAD) are used. For each individual variant, an implied TF is calculated.  Variants with very high implied TF are capped at 200%.  The bandwidth is set to:   
```
bw = min[max(0.2*rawTF,5/wAD, 0.125 * nth highest TF), 3*rawTF,0.1] 
Where n = max (4,8th percentile of variants) 
```
The final TF is set to the highest TF peak with weight > max(3,3%) variants  
#### LOW_COUNT
If the wAD < 50,  and there are at least max(4,0.5% of variants) with 2 or more fragment support, then the ratio of site count with 1 fragment support (N[1]) compared to the site count with >1 fragment support (N[>1]) is used to estimate the support.  Only variants with VCN <= 2* median VCN and sampleDP <=2*wAD are used.   Specifically, the model assumes that a proportion of variants has been dropped from the sample (dropout%) and tests a matrix of combinations of {sampleVAF, dropout%} to see which best explains the ratio N[>1]/N[1].   The combination which best matches the observed ratio is selected.  

^Note that, wAD is weighted for CN and depth in the tissue sample and is defined as  
```
wAD = Σ(i=1->n)[(DPi_cfDNA)2/CNn] / Σ(i=1->n)[DPi_cfDNA/CNi] 
```
The adjusted VAF is used to recalculate the implied tumor fraction estimate for the sample.  If the conditions for neither model are met, then the VAF is not adjusted.    


Finally, the following fields are also calculated 
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


## Output Files
<TO DO>

## Known Issues / Future improvements

<TO DO>

 ## Version History and Download Links


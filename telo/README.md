# TEAL

## Running TEAL

TO DO

### INPUTS

The main inputs for TEAL are a reference and tumor bam.   The following inputs are also required and may be sourced from HMF pipeline outputs (PURPLE and WGS metrcis), or if not specified must be input as values:

Input | Source File | Samples | Default (if not parsed) 
--|--|--|--
Purity | .purple.purity.tsv | Tumor  | 2
Ploidy | .purple.purity.tsv | Tumor | 1
DuplicatePercent | .wgsmetrics | Tumor, Ref | 0
MeanReadsPerKb | .cobalt.gc.median.tsv | Tumor, Ref | Required
GC50ReadsPerKb | .cobalt.gc.median.tsv | Tumor, Ref | (=MeanReadsPerKb)

## Algorithm

There are 4 steps in the algorithm:

### 1. Extraction of telomeric fragments

The bam / cram is sliced for any fragments (mapped or unmapped) where at least 1 read contains 2 or more adjacent canonical telomeric repeats, including supplementary and duplicate reads. 

The output of this step is a Telomere BAM (typically 10,000x smaller than the original BAM) which contains all candidate telomeric fragments. 

### 2. Telomeric annotation 

In this step, each read is annotated according to it’s telomeric and other characteristics. First, any PolyG tails (a common sequencing error which may be confused with G rich telomeric content) are identified (at least 5 consecutive G at end of read) and marked. Then the counts of canonical G orientation (TTAGGG) and C orientation (TAACCC) telomeric repeats are determined and the read is determined to be either ‘C rich’ or ‘G rich’. Finally the counts of other non-canonical telomeric repeats are counted in the orientation determined

### 3. Estimate total telomeric content and length

We must first estimate the total amount of telomeric content (in bases) in the BAM. To do this we count the number of fragments with both reads telomeric and with only one read telomeric. A read is classified as telomeric if it has at least 4 canonical T-type repeats and 5 consecutive telomeric repeats overall. 

Where only 1 read is telomeric the orientation of the telomeric read is important as only C-rich fragments are candidate telomeres, whereas the G-rich) likely represent one end of an interstitial telomeric repeat. As described in TelomereCat (https://www.nature.com/articles/s41598-017-14403-y), we expect interstitial telomeric repeats to be symmetric and have equal numbers of G and C rich reads. Hence we can use the G-rich count to estimate the proportion of the c-rich single read telomeric 
```
Total Telomeric Reads = 2 * Both Telomeric Fragment Count + Single Read Telomeric Fragment C-rich count - Single Read Telomeric Fragment G-rich count
```
To calculate the average telomere length we need to normalise this to the coverage of the genome as a whole. The formula used to normalise is:
```
Mean Telomere Length = Total Telomeric Reads * (1-duplicatePercent) * 1000 / MeanReadsPerKb / GCBiasAdj / 46 
```
The 1000 constant is the COBALT read count window size in bases  and 46 is the number of telomeres in 1 copy of the genome (2 per chromosome).

The GC bias adjustment is set based on the empirical observation of germline telomere content for samples with very high positive or negative GC bias. In practice we calculate GCBias as GC50ReadsPerKb / MeanReadsPerKb and observe that very low values (ie. negative bias) and high values both tend to fit to longer telomere lengths. We apply the following adjustment separately to both tumor and germline samples:

GC50Bias | GCBiasAdj
--|--
0.6 | 1.24
0.65 | 1.06
0.7 | 1.03
0.75 | 1.00
0.8 | 0.98
0.85 | 0.97
0.9 | 0.99
0.95 | 0.99
1 | 0.98
1.05 | 1.00
1.1 | 1.05

For tumor samples, we need to recognise that we are observing a mix of reference and tumor. We use the purity and ploidy obtained from PURPLE, we can solve the following equation for the tumor reference mix. 
```
TumorMixLength = [TumorLength*Purity*Ploidy + RefLength*(1-purity)*2] / [Purity*Ploidy +2*(1-purity)]
```
Rearrangement of this formula gives the following expression for tumor length:
```
TumorLength = [TumorMixLength*[Purity*Ploidy +2*(1-purity)] - RefLength*(1-purity)*2] / purity / ploidy
```
Note this assumes that the stromal component of the tumor has the same telomeric length as the reference sample, but these measurements may differ for different cell types

### 4. Identification of telomeric rearrangements

To identify telomeric rearrangements the telomeric bam is searched for reads with soft clipping in both the germline and tumor sample excluding a small set of blacklisted regions where telomeric reads commonly align (49kb in total of the genome).   Any read that contains a soft clip that contains a segment that is at least 90% match with canonical repeats and at least 12 bases long with the aligned part not containing a segment that is at least 80% match with canonical repeats and at least 12 bases long in the same orientation is considered a candidate telomeric rearrangement site. 

Each candidate location is then inspected and the numbers of split reads and discordant pairs that may support a telomeric rearrangement at the site are counted in both the germline and tumor sample. A read is counted as split read support at a candidate location at the site if it has a soft clip at the precise location. Supplementary reads with hard clipping matching the candidate site are also counted. Where multiple candidate sites with the same orientation and telomeric content orientation are within 50 bases of other then all but the strongest split read support are filtered as ‘DUPLICATE’.

A fragment is counted as a ‘discordant pair’ support if the fragment has a improper pair alignment, neither read overlaps the candidate location, one read faces the site and is within 1000 bases on the non soft clipped side of the candidate location and does not contain a telomeric sequence (as defined above) but has a paired read which does contain a telomeric sequence. Where a discordant pair read may be counted towards multiple candidate locations, the support is assigned to the location nearest to the read alignment. 


## Outputs

To Do



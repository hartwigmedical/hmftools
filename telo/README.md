# TEAL

## Running TEAL

TO DO => Add example usage

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

The following soft filters are applied:

Name | Description | Threshold
--|--|--
maxGermlineSupport | Total # of reads supporting  | <min(5,2% of tumor support) | 
minTelomericLength | Length of longest continuous telomeric content in soft clip or paired read  | >=20
minAnchorLength | Length of longest anchored soft clip or paired read supporting the rearrangement | >=50
minSplitReads | TumorSRTelNoDP + TumorSRTelDPTel + TumorSRNotTelDPTel + TumorSRTelDPNotTel | >=3
minDiscordantPairs | TumorDPTelNoSR + TumorSRNotTelDPTel + TumorSRTelDPTel  | >0
maxCohortFreq | Proportion of samples in cohort with at least 1 rearrangement within +/-50 bases  | <TBD
minTumorMAPQ | Minimum mapq support of locally anchored reads | >=300
duplicate | Must not be a duplicate of another nearby breakend  | =FALSE


The following regions are blacklisted from calling telomeric rearrangements as they are frequently found to have telomeric sequences in recurrent artefacts across samples (currently for hg19/grch37 assembly only)  

Chromosome | Start Position | End Position
--|--|--
1 | 9000 | 11000
1 | 121483000 | 121486000
1 | 249238000 | 249241000
2 | 243152000 | 243154000
3 | 197900000 | 197902000
4 | 9000 | 11000
4 | 191039000 | 191045000
5 | 9000 | 13000
7 | 9000 | 11000
7 | 105741000 | 105743000
8 | 43092000 | 43098000
9 | 9000 | 11000
10 | 42356000 | 42362000
10 | 42377000 | 42401000
10 | 42527000 | 42531000
10 | 42596000 | 42601000
10 | 135523000 | 135526000
12 | 94000 | 97000
12 | 133841000 | 133843000
15 | 102520000 | 102522000
18 | 9000 | 12000
19 | 27731000 | 27739000
19 | 59118000 | 59120000
20 | 62917000 | 62919000
21 | 9704000 | 9706000
21 | 48119000 | 48121000
X | 155259000 | 155262000
MT | 11000 | 14000


## Outputs

The outputs of TEAL is a 'telbam' file (ie a bam restricted to fragments where at least 1 read contains telomeric content), a file which details the estimated  telomeric length and content and finally a file which predicts the telomeric reararrangements
  
### Telomeric Length and content
Column | Description
--|--
SampleId | Unique Id of sample
Type | ‘Tumor’ or ‘Ref’
Purity | From PURPLE (=1 for ref)
Ploidy | From PURPLE (=2 for ref)
GC50ReadsPerKb | From COBALT
MeanReadsPerKb | From COBALT
DuplicatePercent | Estimated proportion of duplicates in file (from WGS metrics)
fullFragments | Count of fragments with both reads classified as telomeric
cRichPartialFragments | Count of fragments with 1 read non telomeric and the other telomeric oriented in a C-rich orientation
gRichPartialFragments | Count of fragments with 1 read non telomeric and the other telomeric oriented in a G-rich orientation
TotalTelomericReads | Count of telomeric reads excluding estimated interstitial repeats
RawTelomereLength | Telomeric reads normalized for total coverage
FinalTelomereLength | Final telomere length adjusted for tumor purity (=Raw length for ref sample)
T-TypeProportion | T-Type repeats / total telomeric repeats
C-TypeProportion | C-Type repeats / total telomeric repeats
G-TypeProportion | G-Type repeats / total telomeric repeats
J-TypeProportion | J-Type repeats / total telomeric repeats
OtherTypeProportion | Other Type repeats / total telomeric repeats


### Rearrangements
Field | Description
--|--
Chromosome | Chromosome of breakend
Position | Position of Breakend
Orientation | 1 = breakend on left side; -1 = breakend on right side
COrGRich | ‘C’ or ‘G’
DistanceToTelomere | Distance to nearest reference telomere in nucleotides
MaxTelomericLength | Longest continuous telomeric segment on any fragment supporting the rearrangement
MaxAnchorLength | Longest locally anchored segment on any fragment supporting the rearrangement
Filter | Either ‘PASS’ or one or more of the below filters
TumorSRTelNoDP | Count of fragments in tumor with 1 read with soft clip supporting the breakend and containing telomeric sequence with the paired read locally anchored
TumorDPTelNoSR | Count of fragments in tumor with 1 read locally anchored but not spanning the breakend and the paired read discordant and containing telomeric sequence
TumorSRTelDPTel | Count of fragments in tumor with 1 read with soft clip supporting the breakend and containing telomeric sequence with the paired read discordant and containing telomeric sequence
TumorSRTelDPNotTel | Count of fragments in tumor with 1 read with soft clip supporting the breakend and containing telomeric sequence with the paired read discordant but not containing telomeric sequence (may indicate telomeric insertion)
TumorSRNotTelDPTel | Count of fragments in tumor with 1 read with soft clip supporting the breakend but NOT containing telomeric sequence with the paired read discordant and containing telomeric sequence
TumorMAPQ | Sum of MAPQ of locally anchored reads in fragments supporting the rearrangement in tumor
GermlineSRTelNoDP | Count of fragments in germline with 1 read with soft clip supporting the breakend and containing telomeric sequence with the paired read locally anchored
GermlineDPTelNoSR | Count of fragments in germline with 1 read locally anchored but not spanning the breakend and the paired read discordant and containing telomeric sequence
GermlineSRTelDPTel | Count of fragments in germline with 1 read with soft clip supporting the breakend and containing telomeric sequence with the paired read discordant and containing telomeric sequence
GermlineSRTelDPNotTel | Count of fragments in germline with 1 read with soft clip supporting the breakend and containing telomeric sequence with the paired read discordant but not containing telomeric sequence (may indicate telomeric insertion)
GermlineSRNotTelDPTel | Count of fragments in germline with 1 read with soft clip supporting the breakend but NOT containing telomeric sequence with the paired read discordant and containing telomeric sequence
germlineMAPQ | Sum of MAPQ of locally anchored reads in fragments supporting the rearrangement in germline
CohortFrequency | Annotated from cohortFreq.bed file

## Known issues and future improvements
* A blacklist bed file is currently only provided for hg19 / grch37 assembly.  
* TEAL represents each breakend independently, but ideally should pair up rearrangements which are supported by the same fragment and represent telomeric insertions
* TEAL should determine a consensus sequence for each telomeric rearrangement




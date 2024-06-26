# TEAL

TEAL measures telomere content, and estimates telomeric length based on WGS BAM input and can be run on a germline only, tumor only or tumor-normal pair.

If a tumor-normal pair is provided, TEAL will also call somatic telomeric rearrangements, ie. breakends linking non telomeric regions of the genome to telomeric content.

## Installation

To install, download the latest compiled jar file from the [download links](#version-history-and-download-links).

TEAL requires Java 11+ to be installed.

## Using TEAL with HMF pipeline
We can run TEAL using output files from HMF pipeline tools.

### Paired germline/tumor mode with HMF pipeline

This is the default mode.
Arguments:

| Argument              | Description                                                                                |
|-----------------------|--------------------------------------------------------------------------------------------|
| reference             | Name of the reference sample                                                               |
| reference_bam         | Path to indexed reference BAM or CRAM file                                                 |
| tumor                 | Name of the tumor sample                                                                   |
| tumor_bam             | Path to indexed tumor BAM or CRAM file                                                     |
| output_dir            | Path to the output directory. This directory will be created if it does not already exist. |
| purple                | Path to PURPLE output. This should correspond to the output_dir used in PURPLE             |
| cobalt                | Path to COBALT output. This should correspond to the output_dir used in COBALT             |
| reference_wgs_metrics | Path to the metrics file of the reference BAM file                                         |
| tumor_wgs_metrics     | Path to the metrics file of the tumor BAM file                                             |
| threads (default = 1) | Number of threads to use                                                                   |
| ref_genome (optional) | Path to the reference genome fasta file. Required only when using CRAM files.              |

Example Usage:

```
java -Xmx16G -cp teal.jar com.hartwig.hmftools.teal.TealApplication \
   -reference COLO829R -reference_bam COLO829R.bam \
   -tumor COLO829T -tumor_bam COLO829T.bam \
   -purple /path/to/COLO829/purple \
   -cobalt /path/to/COLO829/cobalt \
   -reference_wgs_metrics /path/to/COLO829/bam_metrics/COLO829R_WGSMetrics.txt \
   -tumor_wgs_metrics /path/to/COLO829/bam_metrics/COLO829T_WGSMetrics.txt \
   -output_dir /path/to/COLO829/teal \
   -threads 28
```

### Germline only mode with HMF pipeline files

| Argument              | Description                                                                                |
|-----------------------|--------------------------------------------------------------------------------------------|
| reference             | Name of the reference sample                                                               |
| reference_bam         | Path to indexed reference BAM or CRAM file                                                 |
| output_dir            | Path to the output directory. This directory will be created if it does not already exist. |
| cobalt                | Path to COBALT output. This should correspond to the output_dir used in COBALT             |
| reference_wgs_metrics | Path to the metrics file of the reference BAM file                                         |
| threads (default = 1) | Number of threads to use                                                                   |
| ref_genome (optional) | Path to the reference genome fasta file. Required only when using CRAM files.              |

Example Usage:

```
java -Xmx16G -cp teal.jar com.hartwig.hmftools.teal.TealApplication \
   -reference COLO829R -reference_bam COLO829R.bam \
   -cobalt /path/to/COLO829/cobalt \
   -reference_wgs_metrics /path/to/COLO829/bam_metrics/COLO829R_WGSMetrics.txt \
   -output_dir /path/to/COLO829/teal \
   -threads 28
```

## Using TEAL standalone (without HMF pipeline)
If HMF pipeline is not available, then TEAL can be used in standalone mode. Inputs that it extracts from
other pipeline tools will need to be explicited provided.

### Paired germline/tumor mode standalone

| Argument                       | Default                   | Description                                                                                |
|--------------------------------|---------------------------|--------------------------------------------------------------------------------------------|
| reference                      |                           | Name of the reference sample                                                               |
| reference_bam                  |                           | Path to indexed reference BAM or CRAM file                                                 |
| tumor                          |                           | Name of the tumor sample                                                                   |
| tumor_bam                      |                           | Path to indexed tumor BAM or CRAM file                                                     |
| output_dir                     |                           | Path to the output directory. This directory will be created if it does not already exist. |
| reference_duplicate_proportion | 0                         | Proportion of reads that are marked duplicates in the reference sample BAM                 |
| reference_mean_read_depth      |                           | Mean read depth of the reference sample                                                    |
| reference_gc50_read_depth      | reference_mean_read_depth | GC 50 read depth of the reference sample. Defaults to mean read depth if not provided      |
| tumor_purity                   | 1                         | Purity of the tumor sample                                                                 |
| tumor_ploidy                   | 2                         | Ploidy of the tumor                                                                        |
| tumor_duplicate_proportion     | 0                         | Proportion of reads that are marked duplicates in the tumor sample BAM                     |
| tumor_mean_read_depth          |                           | Mean read depth of the tumor sample                                                        |
| tumor_gc50_read_depth          | tumor_mean_read_depth     | GC 50 read depth. Defaults to mean read depth if not provided                              |
| threads                        | 1                         | Number of threads to use                                                                   |
| ref_genome                     |                           | Path to the reference genome fasta file. Required only when using CRAM files.              |

Example Usage:
```
java -Xmx16G -cp teal.jar com.hartwig.hmftools.teal.TealApplication \
   -reference COLO829R -reference_bam COLO829R.bam \
   -tumor COLO829T -tumor_bam COLO829T.bam \
   -reference_duplicate_proportion 0.2284 \
   -reference_gc50_read_depth 21.4 \
   -reference_mean_read_depth 27.1 \
   -tumor_purity 0.6 \
   -tumor_ploidy 1.98 \
   -tumor_duplicate_proportion 0.2766 \
   -tumor_gc50_read_depth 55.5 \
   -tumor_mean_read_depth 57.5 \
   -output_dir /path/to/COLO829/teal \
   -threads 28
```

### Running germline only mode standalone

| Argument                       | Default                     | Description                                                                                |
|--------------------------------|-----------------------------|--------------------------------------------------------------------------------------------|
| reference                      |                             | Name of the reference sample                                                               |
| reference_bam                  |                             | Path to indexed reference BAM or CRAM file                                                 |
| output_dir                     |                             | Path to the output directory. This directory will be created if it does not already exist. |
| reference_duplicate_proportion | 0                           | Proportion of reads that are marked duplicates in the reference sample BAM                 |
| reference_mean_read_depth      |                             | Mean read depth of the reference sample                                                    |
| reference_gc50_read_depth      | reference_mean_read_depth   | GC 50 read depth of the reference sample. Defaults to mean read depth if not provided    |
| threads                        | 1                           | Number of threads to use                                                                   |
| ref_genome                     |                             | Path to the reference genome fasta file. Required only when using CRAM files.              |

Example Usage:

```
java -Xmx16G -cp teal.jar com.hartwig.hmftools.teal.TealApplication \
   -reference COLO829R -reference_bam COLO829R.bam \
   -reference_duplicate_proportion 0.2284 \
   -reference_gc50_read_depth 21.4 \
   -reference_mean_read_depth 27.1 \
   -output_dir /path/to/COLO829/teal \
   -threads 28
```
 
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

$$ Total Telomeric Reads = 2 \times Both Telomeric Fragment Count + Single Read Telomeric Fragment C-rich count - Single Read Telomeric Fragment G-rich count $$

To calculate the average telomere length we need to normalise this to the coverage of the genome as a whole. The formula used to normalise is:

$$ Mean Telomere Length = \frac{Total Telomeric Reads \times (1-duplicatePercent) \times MeanReadLength }{ MeanReadDepth \times GCBiasAdj \times 46 } $$

The 1000 constant is the COBALT read count window size in bases  and 46 is the number of telomeres in 1 copy of the genome (2 per chromosome).

The GC bias adjustment is set based on the empirical observation of germline telomere content for samples with very high positive or negative GC bias.
In practice we calculate GCBias as GC50ReadDepth / MeanReadDepth and observe that very low values (ie. negative bias) and high values both tend
to fit to longer telomere lengths. We apply the following adjustment separately to both tumor and germline samples:

| GC50Bias | GCBiasAdj |
|----------|-----------|
| 0.6      | 1.24      |
| 0.65     | 1.06      |
| 0.7      | 1.03      |
| 0.75     | 1.00      |
| 0.8      | 0.98      |
| 0.85     | 0.97      |
| 0.9      | 0.99      |
| 0.95     | 0.99      |
| 1        | 0.98      |
| 1.05     | 1.00      |
| 1.1      | 1.05      |

For tumor samples, we need to recognise that we are observing a mix of reference and tumor. We use the purity and ploidy obtained from PURPLE, we can solve the following equation for the tumor reference mix. 

$$ TumorMixLength = \frac{ TumorLength \times Purity \times Ploidy + RefLength \times (1-purity) \times 2 }{ Purity \times Ploidy + 2 \times (1-purity) } $$

Rearrangement of this formula gives the following expression for tumor length:

$$ TumorLength = \frac{ TumorMixLength \times [Purity \times Ploidy +2 \times (1-purity)] - RefLength \times (1-purity) \times 2 }{ purity \times ploidy } $$

Note this assumes that the stromal component of the tumor has the same telomeric length as the reference sample, but these measurements may differ for different cell types

### 4. Identification of telomeric rearrangements

To identify telomeric rearrangements the telomeric bam is searched for reads with soft clipping in both the germline and tumor sample excluding a small set of blacklisted regions where telomeric reads commonly align (49kb in total of the genome).   Any read that contains a soft clip that contains a segment that is at least 90% match with canonical repeats and at least 12 bases long with the aligned part not containing a segment that is at least 80% match with canonical repeats and at least 12 bases long in the same orientation is considered a candidate telomeric rearrangement site. 

Each candidate location is then inspected and the numbers of split reads and discordant pairs that may support a telomeric rearrangement at the site are counted in both the germline and tumor sample. A read is counted as split read support at a candidate location at the site if it has a soft clip at the precise location. Supplementary reads with hard clipping matching the candidate site are also counted. Where multiple candidate sites with the same orientation and telomeric content orientation are within 50 bases of other then all but the strongest split read support are filtered as ‘DUPLICATE’.

A fragment is counted as a ‘discordant pair’ support if the fragment has a improper pair alignment, neither read overlaps the candidate location, one read faces the site and is within 1000 bases on the non soft clipped side of the candidate location and does not contain a telomeric sequence (as defined above) but has a paired read which does contain a telomeric sequence. Where a discordant pair read may be counted towards multiple candidate locations, the support is assigned to the location nearest to the read alignment. 

The following soft filters are applied:

| Name               | Description                                                                      | Threshold                   |
|--------------------|----------------------------------------------------------------------------------|-----------------------------|
| maxGermlineSupport | Total # of reads supporting                                                      | <min(5,2% of tumor support) |
| minTelomericLength | Length of longest continuous telomeric content in soft clip or paired read       | >=20                        |
| minAnchorLength    | Length of longest anchored soft clip or paired read supporting the rearrangement | >=50                        |
| minSplitReads      | TumorSRTelNoDP + TumorSRTelDPTel + TumorSRNotTelDPTel + TumorSRTelDPNotTel       | >=3                         |
| minDiscordantPairs | TumorDPTelNoSR + TumorSRNotTelDPTel + TumorSRTelDPTel                            | >0                          |
| maxCohortFreq      | Proportion of samples in cohort with at least 1 rearrangement within +/-50 bases | <TBD                        |
| minTumorMAPQ       | Minimum mapq support of locally anchored reads                                   | >=300                       |
| duplicate          | Must not be a duplicate of another nearby breakend                               | =FALSE                      |


The following regions are blacklisted from calling telomeric rearrangements as they are frequently found to have telomeric sequences in
recurrent artefacts across samples (currently for hg19/grch37 assembly only)  

| Chromosome | Start Position | End Position |
|------------|----------------|--------------|
| 1          | 9000           | 11000        |
| 1          | 121483000      | 121486000    |
| 1          | 249238000      | 249241000    |
| 2          | 243152000      | 243154000    |
| 3          | 197900000      | 197902000    |
| 4          | 9000           | 11000        |
| 4          | 191039000      | 191045000    |
| 5          | 9000           | 13000        |
| 7          | 9000           | 11000        |
| 7          | 105741000      | 105743000    | 
| 8          | 43092000       | 43098000     |
| 9          | 9000           | 11000        |
| 10         | 42356000       | 42362000     |
| 10         | 42377000       | 42401000     |
| 10         | 42527000       | 42531000     |
| 10         | 42596000       | 42601000     |
| 10         | 135523000      | 135526000    |
| 12         | 92000          | 97000        |
| 12         | 133841000      | 133843000    |
| 15         | 102520000      | 102522000    |
| 18         | 9000           | 12000        |
| 19         | 27731000       | 27739000     |
| 19         | 59118000       | 59120000     |
| 20         | 62917000       | 62919000     |
| 21         | 9704000        | 9706000      | 
| 21         | 48119000       | 48121000     |
| X          | 155259000      | 155262000    |
| MT         | 11000          | 14000        |


## Outputs

The outputs of TEAL is a 'telbam' file (ie a bam restricted to fragments where at least 1 read contains telomeric content), a file which
details the estimated  telomeric length and content and finally a file which predicts the telomeric reararrangements
  
### Telomeric Length and content
| Column                | Description                                                                                           |
|-----------------------|-------------------------------------------------------------------------------------------------------|
| purity                | From PURPLE (=1 for ref)                                                                              |
| ploidy                | From PURPLE (=2 for ref)                                                                              |
| gc50ReadDepth         | From COBALT                                                                                           |
| meanReadDepth         | From COBALT                                                                                           |
| duplicateProportion   | Estimated proportion of duplicates in file (from WGS metrics)                                         |
| fullFragments         | Count of fragments with both reads classified as telomeric                                            |
| cRichPartialFragments | Count of fragments with 1 read non telomeric and the other telomeric oriented in a C-rich orientation |
| gRichPartialFragments | Count of fragments with 1 read non telomeric and the other telomeric oriented in a G-rich orientation |
| sampleMixLength       | Telomeric reads normalized for total coverage                                                         |
| telomereLength        | Final telomere length adjusted for tumor purity (=Raw length for ref sample)                          |

### Rearrangements
| Field                 | Description                                                                                                                                                                                                                  |
|-----------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| chromosome            | Chromosome of breakend                                                                                                                                                                                                       |
| position              | Position of Breakend                                                                                                                                                                                                         |
| orientation           | 1 = breakend on left side; -1 = breakend on right side                                                                                                                                                                       |
| cOrGRich              | ‘C’ or ‘G’                                                                                                                                                                                                                   |
| distanceToTelomere    | Distance to nearest reference telomere in nucleotides                                                                                                                                                                        |
| maxTelomericLength    | Longest continuous telomeric segment on any fragment supporting the rearrangement                                                                                                                                            |
| maxAnchorLength       | Longest locally anchored segment on any fragment supporting the rearrangement                                                                                                                                                |
| filter                | Either ‘PASS’ or one or more of the specified filters separated by semi colon                                                                                                                                                |
| inTumor               | true if in tumor                                                                                                                                                                                                             |
| inGermline            | true if in germline                                                                                                                                                                                                          |
| tumorSRTelNoDP        | Count of fragments in tumor with 1 read with soft clip supporting the breakend and containing telomeric sequence with the paired read locally anchored                                                                       |
| tumorDPTelNoSR        | Count of fragments in tumor with 1 read locally anchored but not spanning the breakend and the paired read discordant and containing telomeric sequence                                                                      |
| tumorSRTelDPTel       | Count of fragments in tumor with 1 read with soft clip supporting the breakend and containing telomeric sequence with the paired read discordant and containing telomeric sequence                                           |
| tumorSRTelDPNotTel    | Count of fragments in tumor with 1 read with soft clip supporting the breakend and containing telomeric sequence with the paired read discordant but not containing telomeric sequence (may indicate telomeric insertion)    |
| tumorSRNotTelDPTel    | Count of fragments in tumor with 1 read with soft clip supporting the breakend but NOT containing telomeric sequence with the paired read discordant and containing telomeric sequence                                       |
| tumorMAPQ             | Sum of MAPQ of locally anchored reads in fragments supporting the rearrangement in tumor                                                                                                                                     |
| germlineSRTelNoDP     | Count of fragments in germline with 1 read with soft clip supporting the breakend and containing telomeric sequence with the paired read locally anchored                                                                    |
| germlineDPTelNoSR     | Count of fragments in germline with 1 read locally anchored but not spanning the breakend and the paired read discordant and containing telomeric sequence                                                                   |
| germlineSRTelDPTel    | Count of fragments in germline with 1 read with soft clip supporting the breakend and containing telomeric sequence with the paired read discordant and containing telomeric sequence                                        |
| germlineSRTelDPNotTel | Count of fragments in germline with 1 read with soft clip supporting the breakend and containing telomeric sequence with the paired read discordant but not containing telomeric sequence (may indicate telomeric insertion) |
| germlineSRNotTelDPTel | Count of fragments in germline with 1 read with soft clip supporting the breakend but NOT containing telomeric sequence with the paired read discordant and containing telomeric sequence                                    |
| germlineMAPQ          | Sum of MAPQ of locally anchored reads in fragments supporting the rearrangement in germline                                                                                                                                  |
| cohortFrequency       | Annotated from cohortFreq.bed file                                                                                                                                                                                           |
| 
## Known issues and future improvements
* A blacklist bed file is currently only provided for hg19 / grch37 assembly.  
* TEAL represents each breakend independently, but ideally should pair up rearrangements which are supported by the same fragment and represent telomeric insertions
* TEAL should determine a consensus sequence for each telomeric rearrangement
* TEAL could aslo count relative amount T-Type, C-Type, G-Type and J-Type content per sample (relevant for ALT pathway identification)

# Version History and Download Links
- [1.3.0](https://github.com/hartwigmedical/hmftools/releases/tag/teal-v1.3.0)
  - Ignore consensus reads
- [1.2.2](https://github.com/hartwigmedical/hmftools/releases/tag/teal-v1.2.2)
  - Fix crash when processing none standard chromosome names 
  - Fix crash when processing unpaired reads 
- [1.2.1](https://github.com/hartwigmedical/hmftools/releases/tag/teal-v1.2.1)
  - Fix backward compatibility issue with cobalt v1.16
- [1.2.0](https://github.com/hartwigmedical/hmftools/releases/tag/teal-v1.2.0)
  - Update to match cobalt v1.16+.
  - Use mean read depth and gc50 read depth instead of read count per 1000 bases window.
  - removed reference / tumor mean_reads_per_kb and gc50_reads_per_kb command line arguments in standalone mode
  - added reference / tumor mean_read_depth and gc50_read_depth command line arguments in standalone mode
- [1.1.0](https://github.com/hartwigmedical/hmftools/releases/tag/teal-v1.1.0)
  - Update to use Smith-Waterman instead.
  - Update to parse newer inputs from purple and cobalt
- [1.0.1](https://github.com/hartwigmedical/hmftools/releases/tag/teal-v1.0.1)
  - Fix crash bug in the writing of output file.
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/teal-v1.0)
  - Downgrade error to warning if a read group cannot be completed.
  - Fix loading of cobalt file.
- [1.0_beta](https://github.com/hartwigmedical/hmftools/releases/tag/teal-v1.0_beta)
    - First release

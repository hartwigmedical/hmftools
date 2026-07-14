# PADDLE (Passenger and Driver Dnds Likelihood)

Uses dNdS values to calculate expected driver and passenger mutation rates.

There are 4 steps required to generate dNdS values
1. [Define the cohort to run on](#define-the-cohort-to-run-on)
2. [Gather variants of cohort](#gather-exonic-variants-for-all-samples-in-the-cohort)
3. [Run dNdScv](#run-dndscv)
4. [Generate driver catalog values](#generate-driver-catalog-values)

## Define the cohort to run on 

The first step is to select a set of samples to run dNdS on. The cohort should contain one sample per patient, with the selected sample 
being the one with the highest purity that passes QC.
 
The sample IDs of the cohort are listed in a TSV file, one sample per line:

```tsv
SampleId
SAMPLE_1
SAMPLE_2
```

## Gather exonic variants for all samples in the cohort

The second step is to:
- Gather the exonic somatic variants for each of the samples in the cohort
- Summarise these variants to get the mutational load for each sample

This is done by running `DndsDataBuilder`, providing a directory containing `<SampleId>.purple.somatic.vcf.gz` files:

```bash
java -cp paddle.jar com.hartwig.hmftools.dnds.builder.DndsDataBuilder \
    -sample_id_file sample_ids.txt \
    -output_dir output/ \
    -purple_dir purple/ \
    -threads 8
```

or from the Hartwig SQL database:

```bash
java -cp paddle.jar com.hartwig.hmftools.dnds.builder.DndsDataBuilder \
    -sample_id_file sample_ids.txt \
    -output_dir output/ \
    -db_user ${user} -db_pass ${pass} -db_url ${url}
```

This produces a `dnds_cohort_mut_load.tsv` file with the mutational load counts per sample:

```tsv
SampleId  SnvBiallelic  SnvNonBiallelic  IndelBiallelic  IndelNonBiallelic
SAMPLE_1          8363	          31285             104                612
```

And a `dnds_cohort_variants.tsv.gz` with containing variant data for all samples:

```tsv
SampleId  Chromosome   Position  Ref  Alt  Gene  Biallelic  Hotspot  WorstCodingEffect  CanonicalCodingEffect  RepeatCount
SAMPLE_1        chr7  140753336    A    T  BRAF      false     true           MISSENSE               MISSENSE            0
```

## Run dNdScv

The next step is to run dNdScv on the data collected by step 2. This requires a modified version of the [original dndscv tool](https://github.com/im3sanger/dndscv).
In addition, this step requires a custom HmfRefCDS.RData. Both the modified tool as well as the HMF ref data are kept in Hartwig's internal GCP environment.

To actually run the modified version of dndscv, use the dnds.R script provided alongside the paddle jar,
and pass a working dir that is the directory containing the output of step 2.

The relevant outputs of this step are:
 - DndsMutations.tsv
 - HmfRefCDSCv.tsv

## Generate driver catalog values

Final step is to take the dNdS values from previous step and transform them into something useable by the driver catalog. 

This is done by running the PaddleDriverApplicationKt application:
```
java -cp paddle.jar com.hartwig.hmftools.dnds.calcs.DndsFileBuilder \
    -work_dir /path/to/work_dir
```  

The outputs of this step are:
 - DndsDriverLikelihoodOnco.tsv
 - DndsDriverLikelihoodTsg.tsv
  
These files are ingested into hmf-common resources (hmf-common/src/main/resources/dnds) and are used to calculate driver likelihood of individual mutations of individual samples.

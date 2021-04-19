# PADDLE (Passenger and Driver Dnds Likelihood)

Uses dNdS values to calculate expected driver and passenger mutation rates.

There are 4 steps required to generate dNdS values
1. [Define the cohort to run on](#define-the-cohort-to-run-on)
2. [Gather variants of cohort](#collecting-exonic-variants-for-all-samples-in-the-cohort)
3. [Run dNdScv](#run-dndscv)
4. Update driver catalog values

## Define the cohort to run on 

First step is to create the cohort to run dNdS on. The cohort should not contain more than one sample for any specific patient.
To update the driver catalog values used by Hartwig, the cohort needs to be defined as the Highest Purity Cohort which is the 
best sample for every patient that passes QC. Any smaller cohort could be picked for more specific analyses. 
 
The output from this step should be a tab delimited file as follows:

```
sampleId	purity
SAMPLE090001T	0.48
SAMPLE090002T	0.82
SAMPLE090003T	0.55
SAMPLE090004T	0.42
```

## Collecting exonic variants for all samples in the cohort

Step 2 is to download the exonic somatic variants for each of the samples in the cohort and summarise these variants to get the mutational load for each sample.

This is done by running the PaddleExonicVariantsApplicationKt application:

```
java -cp paddle.jar com.hartwig.hmftools.paddle.PaddleExonicVariantsApplicationKt \
    -output_dir /path/to/outputDir \
    -cohort_tsv /path/to/cohort_tsv_generated_by_step1.tsv \
    -db_user ${user} -db_pass ${pass} -db_url ${url}
```  

This step assumes all samples are present in an HMF patient database that the application can connect to with the credentials provided.

The outputs are:
 - mutationalLoad.tsv: One line per sample describing various mutational loads
 - somatics: A directory containing one file per sample holding all somatic mutations in exonic domain

## Run dNdScv

The next step is to run dNdScv on the data collected by step 2. This requires a modified version of the dndscv tool defined [here](https://github.com/im3sanger/dndscv).
In addition, this step requires a custom HmfRefCDS.RData. Both are available upon request.

To actually run the modified version of dndscv, use the dnds.R script provided alongside the PADDLE jar,
and pass a working dir that is the directory containing the output of step 2.

The outputs are:
 - allSomatics.RData
 - DndsMutations.tsv
 - HmfRefCDSCv.RData
 - HmfRefCDSCv.tsv

## Update driver catalog values

Final step is to take the dNdS values and transform them into something useable by the driver catalog. 

This is done by running the PaddleDndsApplicationKt application:
```
java -cp paddle.jar com.hartwig.hmftools.paddle.PaddleDndsApplicationKt \
    -work_dir /path/to/work_dir
```  

The outputs of this step are:
 - DndsDriverLikelihoodOnco.tsv
 - DndsDriverLikelihoodTsg.tsv
  
These files are ingested into hmf-common resources and are used to calculate driver likelihood of individual mutations of individual samples.

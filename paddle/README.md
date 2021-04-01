# PADDLE (Passenger and Driver Dnds Likelihood)

Uses dNdS values to calculate expected driver and passenger mutation rates.

## Highest Purity Cohort

First step is to create the highest purity cohort. 
This is a list of no more than one sample per patient favouring highest purity when multiple samples are present.
Only samples that have passed QC with a detectable tumor are included. 
By default, there is a minimum purity of 0.2 but this can be overwritten.

```
java -cp paddle.jar com.hartwig.hmftools.paddle.PaddleCohortApplicationKt \
    -db_user username
    -db_pass password
    -db_url mysql://localhost:3307/hmfpatients?serverTimezone=UTC
    -out ~/hmf/analysis/paddle/highestPurityCohort.tsv
    -min_purity 0.2
```  

This will produce a very simple tab delimited file as follows:

```
sampleId	purity
SAMPLE090001T	0.48
SAMPLE090002T	0.82
SAMPLE090003T	0.55
SAMPLE090004T	0.42
```
You can trivially change the cohort to only be for a specific cancer by adding the primary tumor location to the call in DatabaseAccess.readSamplePurityPassingQC().

## HPC Exonic Variants


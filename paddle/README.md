# PADDLE (Passenger and Driver Dnds Likelihood)

Uses dNdS values to calculate expected driver and passenger mutation rates.

There are 4 steps required to refresh dNdS values
1. Create a cohort
2. Gather variants of cohort
3. Run dNdS
4. Update driver catalog values

While the first step is fairly polished the remaining 3 steps require some manual intervention for setting parameters etc.

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
    -out ~/hmf/analysis/paddle/hpc.tsv
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

If you are trying to run locally against datastore don't forget to create an ssh tunnel with a command like:

```
ssh -L 3307:localhost:3306 jon@hmf-datastore
```

## HPC Exonic Variants

Step 2 is to download the exonic somatic variants for each of the samples in the highest purity cohort and summarise these variants to get the mutational load for each sample.

This is done by running the PaddleExonicVariantsApplicationKt application:

```
java -cp paddle.jar com.hartwig.hmftools.paddle.PaddleExonicVariantsApplicationKt
```  
This application does not currently accept arguments. You must overwrite the following lines in PaddleExonicVariantsApplication.kt:

```kotlin
val somaticsDir = "/Users/jon/hmf/analysis/paddle/somatics/"
val cohortFile = "/Users/jon/hmf/analysis/paddle/hpc.head.tsv"
val cohortMutationalLoadFile = "/Users/jon/hmf/analysis/paddle/mutationalLoad.tsv"
```

As this can be a lengthy process, the application is designed to fail (and resume) gracefully. 
It should just be able to pick up where it left off it it gets interrupted.

As with the cohort, a ssh tunnel might be required to connect to the database.


## Run dNdS

This is where things get a little complicated. 

We run a customised older version of dNds that you will need to download and install. You can check the git log to see what changes I made if you ever need/want to use a newer version of dNdS.
I have put it on datastore for you under /data/common/repos/dndscv or there is a tgz file you can download /data/common/repos/dndscv.tgz

We also use a custom HmfRefCDS.RData that I have uploaded to datastore: /data/common/dbs/dnds/HmfRefCDS.RData 

Do not lose these files. They are not easy to replace.

With all the resources downloaded you need to update the (6-ish) file locations in [dnds.R](./src/main/resources/r/dnds.R) and run it. Good luck.

## Update driver catalog values

Final step is to take the dNdS values and transform them into something useable by the driver catalog. 

This is done by running the PaddleDndsApplicationKt application:
```
java -cp paddle.jar com.hartwig.hmftools.paddle.PaddleDndsApplicationKt
```  

This application does not currently accept arguments. You must overwrite the following lines in PaddleDndsApplication.kt:

```kotlin
val path = "/Users/jon/hmf/analysis/paddle" // "dnds4305" || "dnds5441"
val cohortFile = "${path}/mutationalLoad.tsv"
val dndsCVFile = "${path}/HmfRefCDSCv.tsv"
val mutationsFile = "${path}/DndsMutations.tsv"
val hmfToolsRepo="/Users/jon/hmf/repos/hmftools"
```

Note that the output files are designed to be checked into the hmftools repository.


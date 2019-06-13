# SV LINX

A collection of components for analysing and annotating structural variants. 

The key routines are:
- LINE element annotation
- clustering and chaining of variants
- fusion and disruption detection
- driver-gene annotation

Data is loaded from the HMF patients data and other reference files, and each analysis routine writes out applicable annotation results.

## Contents

* [Configuration](#configuration)
* [Defintions and Key Classes](#defintions-and-key-classes)
* [Pre-Clustering Routines](#pre-clustering-routines)
* [Clustering](#clustering)
* [Chaining](#chaining)
* [Fusion Detection](#fusions-analysis)
* [Driver Gene Analysis](#driver-gene-annotation)
* [Version History](#version-history)

## Dependencies

Required DB tables:
* copyNumber
* structuralVariant
* svAnnotation, svCluster and svLink
* svBreakend and svFusion
* purity - can be empty but needs to exist (but see note about QC below)
* geneCopyNumber - can be empty but needs to exist 
* driverCatalog - can be empty but needs to exist

## Configuration
All values are optional unless otherwise specified.

Example command and arguments:

```
java -jar sv-linx.jar 
    -sample SAMPLE_ID 
    -db_url [db_url] -db_user [username] -db_pass [password] 
    -output_dir /path_to_sample_data/ 
    -fragile_site_file fragile_sites.csv 
    -line_element_file line_elements.csv 
    -replication_origins_file heli_rep_origins.bed 
    -viral_hosts_file viral_host_ref.csv 
    -gene_transcripts_dir /path_to_ensembl_data_cache/ 
    -check_fusions 
    -fusion_pairs_csv knownFusionPairs.csv 
    -promiscuous_five_csv knownPromiscuousFive.csv 
    -promiscuous_three_csv knownPromiscuousThree.csv 
    -write_vis_data 
    -log_debug
```


#### Core
Argument  | Description
---|---
sample  | Either a specific sample, list of samples separated by ';' or blank or '*' to process all samples passing QC from the database
output_dir | Required, directory where all output files are written
database connectivity | db_user, db_pass and db_url

#### Modes and Routines
Argument  | Description
---|---
check_drivers | run driver annotation logic
check_fusions | discover and annotate gene fusions

#### Reference files
Argument  | Description
---|---
fragile_site_file | list of known fragile sites
line_element_file | list of known LINE elements
replication_origins_file | replication timing input
viral_hosts_file | list of known viral hosts
gene_transcripts_dir | directory for Ensembl reference files
sv_data_dir | directory for SV flat-file when not using database
purple_dir | directory for purple and purity data files when not using database

#### Clustering
Argument  | Description
---|---
proximity_distance | (default = 5000), minimum distance to cluster SVs 
chaining_sv_limit | threshold for # SVs in clusters to skip chaining routine (default = 2000)
annotations | list of additional annotations on clusters, links and SVs separated by ';', use 'ALL', (default is none)

#### Driver Gene Annotation
Argument  | Description
---|---
write_gcn_data | extract GeneCopyNumber data from the and write it to file for faster usage on subsequent runs
gcn_data_file | load GeneCopyNumber data from file instead of retrieving from DB (an option where the GCN table is slow to access due to its size)

#### Fusions Analysis
Argument  | Description
---|---
sample_rna_file | RNA data from StarFusion (see section on Fusions for more info)
skip_fusion_output | fusion data is not written to file
fusion_gene_distance | distance upstream of gene to consider a breakend applicable (default = 100K)
restricted_fusion_genes | restrict fusion search to specified genes, separated by ';'
no_fusion_phase_match | find fusions without requiring a up & down stream phase match
log_reportable_fusion | only log reportabl fusions
fusion_pairs_csv, promiscuous_five_csv, promiscuous_three_csv | reference files for known fusions andand promiscuous genes


#### Logging
Argument  | Description
---|---
write_vis_data | write output to for generation of Circos clustering and chaining plots
log_debug | logs in verbose mode
log_cluster_id | log data for a specific cluster
log_sv_id | log data for specific SV


## Defintions and Key Classes
Class or Term | Description
---|---
SV | As sourced from the StructuralVariant database table and annotated with various other information
Breakend | The start or end of an SV defined by its position and orientation
Simple SV | DEL, DUP or INS - ie not a BND, INV or SGL
Cluster | A set of SVs grouped by one or more clustering rules such as proximity. A cluster represents a set of SVs which are considered to have occurred in the genomic event.
Linked Pair | 2 SVs from either 2 breakends facing away, forming a deletion bridge (DB), or facing towards each other, forming a templated insertion (TI) 
Inferred vs Assembled Link | Assembled links are those formed using GRIDSS assembly data. Inferred links are speculative links. Both cannot be shorter than the anchor distance specified per breakend by GRIDSS if known. 
Chain - 1 or more TI link pairs linked to form a continuous chain


## Pre-clustering Routines
Prior to clustering the following annotation routines are completed:
- SVs are marked as having either 1 or both breakends in fragile sites  
- SVs are marked as having either 1 or both breakends in known LINE elements
- breakends within or upstream (by default by 100K beses) of genes are marked with all valid transcripts
- a 'long DEL and DUP' length is calculated for each sample (see below)
- SVs are marked with their replication timing
- LOH events are identified and the linking breakends marked
- a min and max ploidy is calculated from copy number data
- mark each SV breakend's nearest neighbouring breakend (other than its own other breakend), and whether it overlaps it
- mark any duplicate breakend
- mark any 2 breakends which form a deletion bridge regardless of the length, and allowing for overlap up to the distance specified by the anchor distances of the breakends

#### LOH Event Identification
A LOH segment is started by the minor allele copy number dropping below 0.5 and rising above this at the other end.

Skipped when ...


#### Suspect LINE Elements
In addition to the Known LINE elements loaded from the reference file, the following method is used to identify likely or 'Suspect' line elements:
- if 2 or more BNDs are within 5K bases with at least one not forming a DB of < 30 bases AND at least one SV within 5kb having an INS sequence containing at least 11 repeated As or Ts.
- OR at least 1 BND with a remote SGL forming a DB of < 30 bases AND EITHER at least one SV also within 5kb OR the remote SGL HAVING an INS sequence containing at least 11 repeated As or Ts. 

Subsequently mark all clusters as type LINE if:
- it contains a suspected LINE element
- if every variant in the cluster is part of a KNOWN line element AND at LEAST one of the variants is a BND
- If it is a SGL or SGL=2 with at least one PolyA insert sequence
- If it is an INS with at least one PolyA insert sequence

#### Long DEL and DUP length calculation
The following method is used to determine a length threshold for what is then considered a 'long' DEL or DUP:
- find all DUPs and DELs on arms with no inversions
- exclude the 5 longest lengths
- interpolate/normalise this figure to number of simple arms vs max (41, being 46 - 5 with no short arm)

#### Ploidy Recalculation
Theoretically Ploidy equals the CopyNumberChange and both the start and end breakend, but in practice they can differ substantially. 

The following routine calculates a ploidy estimate and uncertainty to aid with downstream comparison of SV ploidies.

##### 1. Estimate an uncertainty for copy number change at each breakend
For this use the principle that the uncertainty in copy number change is driver primarily by the uncertainty in the copy number of the least confident adjacent copy number region which in turn is driven primarily by the number or read depth windows used to estimate the length of the adjacent regions.

Hence we use the following formula to

	CNChangeUncertainty = MAX(maxAdjacentCopyNumber * BaseRelativeUncertainy [0.1],BaseAbsoluteUncertainty [0.15])+  MAX(AdditionalAbsoluteUncertainty [0.4],AdditionalRelativeUncertainty [0.15]*maxAdjacentCopyNumber)/SQRT(minAdjacentDepthWindowCount)

If the mindAdjacentDepthWindowCount = 0, then this means the segment is inferred by the SV ploidy in PURPLE already and no copy number estimate is calculated.

For the special case of foldback inversions, if the flanking depth window counts are both higher than the internal depth window count, a single copy number change observation of half the combined copy number change is made with the confidences determined from the flanking windows

##### 2. Estimate an uncertainty for ploidy

The raw ploidy of the SV is already estimated in PURPLE by multiplying the purity adjusted VAF of the SV by the copyNumber at each breakend.    The VAF estimate depends ultimately on the measured readcount of supporting tumor fragments which is a binomial distribution.

To estimate the uncertainty in the PLOIDY, we estimate the 1% and 99% confidence intervals of the true read count from the  observed read count and then calculate the ploidy uncertainty as half the relative range of the confidence interval.

	Ploidy Uncertainty = Ploidy * (ReadCount99.5%CI-ReadCount0.5%CI) / 2 / ObseverdReadCount


##### 3. Average the 3 ploidy predictions and estimate a consolidated uncertainty

To we weigh the observations by the inverse square of their estimated uncertainties:

	 consloidatedPloidy =  SUM[Observation(i)*(1/Uncertainty(i)^2)] / Sum[1/Uncertainty(i)^2]

The combined uncertainty is estimated as the square root of the weighted sum of squares of the difference between the final ploidy estimate and each individual estimate, but capped at a minimum of half the input uncertainty.   Ie. 

	consolidatedUncertainty = SQRT(countObservations /  (countObervations-1) * SUM[1/Uncertainty(i)^2*(MAX(Observation(i)-consolidatedPloidy,Uncertainty(i)/2))^2] / Sum[1/Uncertainty(i)^2] )


## Clustering 

All SVs within a sample are grouped into clusters by the rules described below.

SVs which are excluded from some or all of the clustering rules are:

Type or Condition | Description
---|---
LowQual | Any variant (excluding INS) with copy number change at both ends < 0.5 unless clustered by proximity
Low VAF SGLs | Any SGL with low VAF (< 0.1)
Poly C/G SGLs | Any SGL with a poly C or G moitf (repeat count 16)
Equiv SGLs | Any SGL marked as 'eqv' in GRIDDS assembly data or a duplicate of another breakend
Overlapping DELs | A DEL will not be clustered with another DEL which overlaps it
LINE | Clusters with breakends marked as LINE elements are excluded from all except proximity-clustering rules

#### Simple Clustering Rules 

##### 1. Proximity
Any SV with a breakend within the specified proximity_distance (defaults to 5K bases) of another SV's breakend causes the SVs to be clustered.  

##### 2. Long DELs, DUPs and INVs
Merge any clusters or SVs where a pair of SVs of type DEL, DUP or INV have lengths exceeding the Long DEL-DUP threshold and overlap each other.

##### 3. Unresolved SGLs
Merge clusters with 1 unresolved SGL with the following rules:
- 2 x cluster-1s with SGLs that are each other's nearest neighbours
- only apply a rule between the 2 closest breakends at the exclusions of the cluster on their other end
- unless the other breakend is a short, simple SV

##### 4. LOH events
The 2 breakends forming an LOH are clustered if not a simple DEL.


#### Chaining of Simple Clusters 
Prior to subequent merging of clusters by more complex rules and interaction of SVs, some basic chaining and cluster identification is performed:

##### Simple Cluster Resolution
Clusters comprised of a between 1-3 simple SVs are marked as 'Simple' and by default excluded from subsequent clustering and chaining.

##### Simple Chaining
Clusters not marked as 'Simple' are put through a simple chaining routine:
1. If the cluster has consistent ploidy, has 3 or less SVs and has consistent (see note) breakends, the consider it 'Simple', otherwise 'Completx'.
2. For simple clusters, form a chain allowing inferred links as well as assembled links. For complex clusters, only form chains from assembled links.
3. For simple clusters, identify synthetic DELs and DUPs from the following pairs of SVs:
    - a DEL-DUP pair, 2 non-overlapping DELs or 2 DUPs
    - 2 BNDs connecting the same chromosomal arms
    - 2 INVs
    - 2 SGLs

NOTE: Consistent here means that the sum of breakends facing a telomere equals those facing a centromere.

##### Synthetic DELs, DUPs and INSs
Conditions to form a synthetic DEL or DUP are as follows:

From 2 SGLs:
1. Facing breakends within the minimum TI length (ie using anchor distance) form a INS
2. Facing breakends outside the minimum TI length form a DUP
3. Facing breakends outside the minimum TI length form a DEL

For the remaining pairs the SVs must form at least 1 TI and it if is inferred it cannot traverse another non-Simple cluster.

From 2 BNDs:
1. If both breakend pairs form DBs, mark the cluster as a Reciprocal Translocation
2. If 2 TIs are formed, consider the assembled or shorter one as 'external'
3. If the 2 breakends not forming the external TI form a DB, the cluster is a DEL with External TI, otherwise it is a DUP with External TI

From 2 INVs or a pair comprised of DELs and DUPs:
1. If 2 DBs are formed, form a DEL with Internal TI
2. If a TI is enclosed by DB, form a DEL with Internal TI
3. If a TI is external to a DB, form a DEL with External TI
4. If a TI is external to facing breakends, form a DUP with External TI
5. If a TI overlaps facing breakends, form a DUP with Internal TI
6. If the pair of SVs forms 2 TI, take the shorter of the 2 as the internal or external TI. 
 

#### Complex Clustering Rules 

##### 5. Foldbacks 
A pair of breakends are marked as forming a foldback if:
- the breakends are consecutive
- a DB is permitted on the backmost of 2 forward facing consecutive breakends, but not on the foremost breakend
- copy number between the breakends is consistent
- the breakends are linked by an assembled chain <= 5K bases
- a BND is assembled on one end to both ends of another variant (a replicated translocation) 

Subsequently merge any 2 clusters with a foldback on the same chromosomal arm.

Additionally merge any cluster with a foldback to another cluster if the other cluster has a non-simple SV breakend which faces the foldback < 5M bases away, with matching ploidy.
 

##### 5. Common Arms
Merge any 2 clusters if they 2 touch the same 2 chromosomal arms with the following conditions:
- SVs in LINE elements are ignored
- SVs which link the 2 arms but are in a short (< 1K) TI are ignored

##### 6. Arm End Ploidy 


##### 7. LOH Breakend Resolution


##### 8. Net CopyNumber Resolution
For the first and last uninterrupted footprint of each cluster on each chromosome, calculate the minimal number of telomeric/centromeric facing ploidy that is required to explain the orientation and ploidy of the breakends limited to the major allele ploidy immediately flanking the cluster.

If this exceeds the telomeric / centromeric major allele ploidy then search for other (non-resolved) clusters that could explain the drop in ploidy. If there is only one cluster that can explain the full change in major allele ploidy cluster with that, else choose the nearest cluster.



## Chaining 



### Post-Clustering Annotations

##### Shattering



##### Double Minutes



##### Link Analysis



## Fusion Detection




## Driver Gene Analysis


## Other SV Analysis routines

Linx has routines for other SV-related analyses:
* statistical co-occurence routines, eg gene to bucket
* analyse multiple biopsy samples for shared vs private overlaps
* run simulations of shattering

#### Copy Number Analysis
Argument  | Description
---|---
run_cn_analysis | required to run the following routines
sv_ploidy_file | load calculated ploidy per SV from file
loh_file | load LOH event data from file
write_ploidy_data | calculate and write calcualated ploidy data to file, with optional config 'verbose_ploidy_data' to show intermediary calc values

#### Multiple Biopsy Analysis
Argument  | Description
---|---
patient_ids_file | linking SampleIds by PatientId
sva_svs_file | SVs data (from SV Lynx output) for each SampleId specified in the patient_ids_file


## Version History

- 1.0
    - initial version with clustering, chaining, driver gene annotation and fusion detection
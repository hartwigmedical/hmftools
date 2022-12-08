# CIDER
Using WTS/WGS or targeted data, we determine a comprehensive list of CDR3 sequences for each of the IG and TCR loci including an abundance estimate.  Notably this includes incomplete rearrangements and IGK-KDE deletions.

The intended purposes of this are the following: 
- For B cell & T cell tumors, the IG / TCR loci can be used as a proxy for the number of clones in the tumor. The CDR3 sequence can also be used for MRD detection.  
- For other tumors, the abundance estimate can be used to estimate immune infiltration and may give insight into diversity, and/or recurrent and cancer specific T-cell clones. 
- We may be able to use epitope-TCR binding prediction tools (such as TCRex, NetTCR, Repitope) to predict binding of specific sequences 
- With very deep targeted sequencing we could determine the full IG/TCR receptor repertoire (eg. TCR beta diversity is estimated at 1-3M distinct sequences per individual). Diversity and evenness are also proposed as important characteristics 

## Usage

```
java -Xmx16G -cp cider.jar com.hartwig.hmftools.cider.CiderApplication \
   -sample COLO829T \
   -bam COLO829T.bam \
   -output_dir /path/to/COLO829/cider \
   -ref_genome_version 37 \
   -write_cider_bam \
   -ensembl_data_dir /path_to_ensembl_data_cache/ \
   -threads 16
```

## Algorithm
### Anchor sequences and coordinates 
To create reference data, we have queried from IMGT (https://www.imgt.org/genedb/) to get all sequences for species Homo Sapiens and (separately) for Molecular Component: IG and TR.   Then we select all query results choosing “F+ORF+in-frame P nucleotide sequences with IMGT gaps”.    Following this, for all genes which exist in ensembl (separately for 37 and 38) we have deterimined reference genome anchor coordinates for each gene. 

Specifically, for each V and each J component define a 30 anchor region and obtain the sequence and genome coordinates for all alleles:

| Gene component | Anchor region ( from the human_IMGT+C.fa, must be at least 25 bases)                                                                                                                                                 | 
|----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| IGHV           | Base 283-312.                                                                                                                                                                                                        |  
| IGHJ           | 30 base sequence starting with TGGGG (J-TRP)                                                                                                                                                                         | 
| IGKV           | Base 283-312.                                                                                                                                                                                                        | 
| IGKJ           | Sequence starting with “TTTG” or “TTCG” starting between 27 and 33.                                                                                                                                                  |  
| IGK-KDE        | Specific anchor coordinates (2:88832743-88832772 in hg38) rom the K Del region are checked so that we can detect IGK deletion type rearrangements                                                                    |
| IGLV           | Base 283-312.                                                                                                                                                                                                        |  
| IGLJ           | Sequence starting with “TTTG” or “TTCG” starting between 27 and 33                                                                                                                                                   | 
| TRAV           | Base 283-312.                                                                                                                                                                                                        |  
| TRAJ           | Sequence starting with “TTTG” or “TTCG” starting between 27 and 33 bases from the end of the fa sequence (​​note that TRAJ35,TRAJ33 & TRAJ55 will be excluded because they don’t appear to have the conserved J-PHE) | 
| TRBV           | Base 283-312.                                                                                                                                                                                                        |  
| TRBJ           | Sequence starting with “TTTG” or “TTCG” starting between 27 and 33                                                                                                                                                   | 
| TRDV           | Base 283-312                                                                                                                                                                                                         | 
| TRDJ           | Sequence starting with “TTTG” or “TTCG” starting between 27 and 33                                                                                                                                                   | 
| TRGV           | Base 283-312                                                                                                                                                                                                         | 
| TRGJ           | Sequence starting with “TTTG” or “TTCG” starting between 27 and 33                                                                                                                                                   |

The last 3 bases of the V anchor and the first 3 bases of the J anchor ,corresponding to the conserved C on the V side and W/F on the J side are used as the anchor coordinates. The D region is generally too short to find aligned reads and is ignored. 

### Bam Extraction 

From the bam retrieve all reads and their mates which overlap the anchor coordinates. A V+J bam is created with the extracted reads. We also obtain reads that overlap the first base of each gene in the constant region with a J anchor blosum62 match in soft clipping region (10 AA), along with any unmapped reads with read pairs mapped within 1kb upstream of V gene  or 1KB downstream of J gene or which map to a constant region and have a blossum62 match to a J or V anchor (10AA) of the same gene locus . 

### VDJ consensus candidate sequences 

Separately for V and J aligned / anchored reads, determine the 30 base anchor sequence + any candidate CDR3 sequence starting from the bases immediately following the conserved V-CYS r J-PHE/J-TRP location. Determine a minimal set of consensus sequences separately for each of the V and J side by collapsing sequences which match (trimmed for bases with qual < 25) into a single consensus sequence. Note that TRA and TRD sequences may also match together. PolyG tails of 6 or more consecutive Gs (and the prior 5 bases) are stripped from the sequence before making the consensus sequence.  Additionally num_trim_bases is set to > 0 then the specified number of bases is always trimmed from every read prior to creating the candidate sequences  

If a sequence can be collapsed to multiple longer sequences, then greedily allocate it to the most supported sequence. The total base qual supporting each base is retained for later matching). 

### Identify anchors and call CDR3 sequences  

For merged sequences, find a V and J anchor simply read out the CDR3 sequence between the V-CYS and J-PHE/J-TRP anchors.  

For each V anchored only consensus sequence we search for candidate J-PHE/J-TRP anchor sequences. To do this we compare each complete 10 amino acid kmer downstream of the V-CYS sequence to the set of known 10 amino acid J anchor sequences by summing the log likelihoods from the BLOSUM62 substitution matrix. Truncated partial anchor sequences of 1 or more amino acids are also checked for the final 1-9 amino acids of the consensus sequence. A similarity score is calculated for each anchor sequence as follows: 

Similarity Score = 3 * Amino Acid length -6 - SUM[Self BLOSUM62 - Anchor BLOSUM62]] 

If the max similarity score to any anchor sequence is greater than 0 we deem it to be a CDR3 sequence. If 2 candidate anchors share the same score, then rank first by inframe and then by whichever CDR3 sequence is closest to 13 amino acids in length. For each J anchored only read/fragment we similarly search for candidate V-CYS anchor sequences in the consensus sequence 

### Collapse consensus sequences 

A sequence is collapsed into another sequence if it is identical at all bases with high quality support across both anchors and candidate CDR3 sequence. Sequences may also be collapsed where they have a single base difference with up to 2 reads of high quality support.  

Each V only anchored read is also checked for partial overlap with each J only anchored read. If they exactly match with more than 15 nucleotides of exact match (allowing for low quality bases) then they are also collapsed into a single sequence. 

### Filter and output VDJ sequences 
Each collapsed sequence is either marked as PASS or one or more of the following filters 
- **NO_V_ANCHOR** - No candidate V anchor found 
- **NO_J_ANCHOR** - No candidate J anchor found 
- **DUPLICATE** - CDR3 nt sequence is identical to another sequence with more support (different anchors) 
- **CDR3_DELETED** - A V and J anchor are found, but the CDR3 portion of the sequence (including conserved C,W,F) is fully deleted
- **MAX_LENGTH** - CDR3 nt sequence must be less than 40 AA in length 
- **MIN_LENGTH** - CDR3 nt sequence must be at least 5 AA in length (including anchor C & W/F)
- **MATCHES_REF** - NonSplitRead+vNonSplitReads >=2 AND either vAlignedReads or jAlignedReads=0.

Note that sequences with "no anchor" may represent partial rearrangements.
 
The full set of fields output are:

| Field                | Explanation                                                                                                                                         | 
|----------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|
| CDR3Seq              | CDR3 nucleotide sequence. If either the V or J anchor is missing only the first 63 bases of sequence are shown                                      | 
| CDR3aa               | CDR3 aa sequence. If either the V or J anchor is missing only the first 63 bases of sequence are shown                                              | 
| Filter               | PASS if viable CDR3 sequence or one or more filter reasons  (see above)                                                                             | 
| minHighQualBaseReads | number of reads in the least supported base in the CDR3 region or for the first 63 bases of the candidate CDR3 sequence if only one anchor is found |
| assignedReads        | Total reads assigned to candidate sequence.                                                                                                         | 
| vAlignedReads        | # of reads initially aligned to V gene                                                                                                              | 
| jAlignedReads        | # of reads initially aligned to J gene                                                                                                              | 
| inFrame              | CDR3 sequence is inframe {T/F}                                                                                                                      | 
| containsStop         | CDR3 contains stop codon {T/F}                                                                                                                      | 
| vType                | {IGKV;IGLV;IGHV;TRAV;TRBV;TRDV;TRGV}                                                                                                                | 
| vAnchorEnd           | Position of V anchor end in nt sequence                                                                                                             | 
| vAnchorNT            | V anchor sequence in nt                                                                                                                             | 
| vAnchorTemplateSeq   | Best scoring V template anchor in nt (or null if read aligned to V anchor)                                                                          | 
| vAnchorAA            | V anchor sequence in AA                                                                                                                             | 
| vAnchorTemplateAA    | Best scoring V template Anchor in AA (or null if read aligned to V anchor)                                                                          | 
| vSimilarityScore     | Blosum similarity score for template anchor (or null if read aligned to V anchor)                                                                   | 
| vNonSplitReads       | Count of reads supporting sequence with at least 30 aligned bases either side of last base of conserved C                                           | 
| jType                | {IGHJ;IGKJ;IGK-KDE,IGLJ;TRAJ;TRBJ;TRDJ;TRGJ}                                                                                                        | 
| jAnchorEnd           | Position of J anchor end in nt sequence                                                                                                             | 
| jAnchorSeq           | J anchor sequence in nt                                                                                                                             | 
| jAnchorTemplateSeq   | Best scoring J template anchor in nt (or null if read aligned to J anchor)                                                                          | 
| jAnchorAA            | J anchor sequence in AA                                                                                                                             | 
| jAnchorTemplateAA    | Best scoring J template Anchor in AA (or null if read aligned to J anchor)                                                                          | 
| jSimilarityScore     | Blosum62 similarity score for template anchor (or null if read aligned to J anchor)                                                                 | 
| vNonSplitReads       | Count of reads supporting sequence with at least 30 aligned bases either side of first base of conserved W/F                                        | 
| vdjSeq               | Full consensus sequence in nucleotides                                                                                                              | 
| support              | Counts of high quality base support at each nucleotide (radix-36 ASCII encoded)                                                                     | 
| SampleType           | "dna" or "rna"                                                                                                                                      | 
| cohortFrequency      | TO DO                                                                                                                                               | 


## Limitations / Future improvements
  
### Bam extraction:
- **Reads mapped to other locations** - We only use reads where the alignment overlaps a known V or J anchor sequence coordinate which means the program is fast. We could also look for more reads with sequences that precisely or partially match known anchors but which have not been mapped to the expected locations.    
- **Mate overlap** - Where fragment lengths are short the reads may overlap (particularly relevant for RNA). For each extracted read pair test for overlap by searching for an exact match for the innermost 10 bases of each read (allowing for differences if base quality < 25). If a match is found then check that the full overlapped region is identical (again allowing for base quality trimming). Create a consensus sequence for the 2 reads, using the highest base quality where the sequences differ.  
- **Fragments with both reads unmapped reads** - these are not queried and extracted. 

### CDR3 calling:
- **Full receptor sequence** - We could assemble outwards from the CDR3 to predict the full receptor sequence.  
- **PON** - We should filter sequences found in a large number of samples 
- **Error tolerance in collapsing** - We collapse sequences with up to 1 high quality sequencing difference across the anchors + CDR3 sequence. We still see a small number of artefacts from very highly supported sequences which could be cleaned up further. 
- **Extension of incomplete TCR** - For TCR regions it may be possible to predict a full CDR3 sequence from the germline using a parital frgament.  For IG this is likely dangerous due to hypermutation 
- **Multiple CDR3s in consensus sequence** - A single consensus sequence may have 2 anchor locations that lead to plausible high scoring CDR3 sequences. Currently we choose the highest scoring, but both could be functional. 


## Post CIDER annotation script using BLAST

Included here is a python annotation script [cider_blastn.py](./src/main/resources/cider_blastn.py) that uses [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK62051/def-item/blast/)
to query each sequence CIDER found against the human genome (GCF_000001405.39_top_level). It uses this information to assign V, D, J alleles and also
weed out false positives.

To run the tool
```
python cider_blastn.py \
    --in_tsv=COLO829.cider.vdj_seq.tsv \
    --out_tsv=COLO829.cider.vdj_seq_blastn.tsv \
    --blastn=/bin/blastn \
    --ensembl=/path/data/ensembl_data_cache/38/ensembl_gene_data.csv
```

This will create a new file `COLO829.cider.vdj_seq_blastn.tsv` with the original file `COLO829.cider.vdj_seq.tsv` plus the following
additional fields added:

| Field                                 | Explanation                                                                    | 
|---------------------------------------|--------------------------------------------------------------------------------|
| vGene, dGene, jGene                   | The V, D or J gene alleles that this sequence is aligned to                    | 
| vPIdent, dPIdent, jPIdent             | The align sequence % identity with the V, D or J gene                          | 
| vAlignStart, dAlignStart, jAlignStart | Start of the alignment with the V, D or J gene                                 | 
| vAlignEnd, dAlignEnd, jAlignEnd       | End of the alignment with the V, D or J gene                                   |
| blastnFullMatch                       | True if this whole sequence matches a part of the genome, False otherwise      | 
| blastnStatus                          | SKIPPED_BLASTN, V_D_J, V_J, V_D, D_J, V_ONLY, D_ONLY, J_ONLY, NO_REARRANGEMENT | 

### setting up BLAST+

Follow the instruction in https://www.ncbi.nlm.nih.gov/books/NBK1762/ . And set up the `human_genome` blast DB:

```
$ perl $BLAST_INSTALL/bin/update_blastdb.pl --passive --decompress human_genome
```

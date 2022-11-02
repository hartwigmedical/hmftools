# CIDER
Using WTS/WGS or targeted data, we determine a comprehensive list of CDR3 sequences for each of the IG and TCR loci including an abundance estimate. 

The intended purposes of this are the following: 
- For B cell & T cell tumors, the IG / TCR loci can be used as a proxy for the number of clones in the tumor. The CDR3 sequence can also be used for MRD detection.  
- For other tumors, the abundance estimate can be used to estimate immune infiltration and may give insight into diversity, and/or recurrent and cancer specific T-cell clones. 
- We may be able to use epitope-TCR binding prediction tools (such as TCRex, NetTCR, Repitope) to predict binding of specific sequences 
- With very deep targeted sequencing we could determine the full IG/TCR receptor repertoire (eg. TCR beta diversity is estimated at 1-3M distinct sequences per individual). Diversity and evenness are also proposed as important characteristics 

## Usage

<TO DO: Hong>

## Algorithm
### Anchor sequences and coordinates 
To create reference data, we have queried from IMGT (https://www.imgt.org/genedb/) to get all sequences for species Homo Sapiens and (separately) for Molecular Component: IG and TR.   Then we select all query results choosing “F+ORF+in-frame P nucleotide sequences with IMGT gaps”.    Following this, for all genes which exist in ensembl (separately for 37 and 38) we have deterimined reference genome anchor coordinates for each gene. 

Specifically, for each V and each J component define a 30 anchor region and obtain the sequence and genome coordinates for all alleles:
  
Gene component | Anchor region ( from the human_IMGT+C.fa, must be at least 25 bases) 
--|--
IGHV | Base 283-312.  
IGHJ | 30 base sequence starting with TGGGG (J-TRP) 
IGKV | Base 283-312. 
IGKJ | Sequence starting with “TTTG” or “TTCG” starting between 27 and 33 
IGLV | Base 283-312.  
IGLJ | Sequence starting with “TTTG” or “TTCG” starting between 27 and 33 
TRAV | Base 283-312.  
TRAJ | Sequence starting with “TTTG” or “TTCG” starting between 27 and 33 bases from the end of the fa sequence (​​note that TRAJ35,TRAJ33 & TRAJ55 will be excluded because they don’t appear to have the conserved J-PHE) 
TRBV | Base 283-312.  
TRBJ | Sequence starting with “TTTG” or “TTCG” starting between 27 and 33 
TRDV | Base 283-312 
TRDJ | Sequence starting with “TTTG” or “TTCG” starting between 27 and 33 
TRGV | Base 283-312 
TRGJ | Sequence starting with “TTTG” or “TTCG” starting between 27 and 33 
  
The last 3 bases of the V anchor and the first 3 bases of the J anchor ,corresponding to the conserved C on the V side and W/F on the J side are used as the anchor coordinates. The D region is generally too short to find aligned reads and is ignored. 

### Bam Extraction 

From the bam retrieve all reads and their mates which overlap the anchor coordinates. A V+J bam is created with the extracted reads. We also obtain reads that overlap the first base of each gene in the constant region with a J anchor blosum62 match in soft clipping region (10 AA), along with any unmapped reads with read pairs mapped within 1kb upstream of V gene  or 1KB downstream of J gene or which map to a constant region and have a blossum62 match to a J or V anchor (10AA) of the same gene locus . 

### VDJ consensus candidate sequences 

Separately for V and J aligned / anchored reads, determine the 30 base anchor sequence + any candidate CDR3 sequence starting from the bases immediately following the conserved V-CYS r J-PHE/J-TRP location. Determine a minimal set of consensus sequences separately for each of the V and J side by collapsing sequences which match (trimmed for bases with qual < 25) into a single consensus sequence. PolyG tails of 6 or more consecutive Gs (and the prior 5 bases) are stripped from the sequence before making the consensus sequence.  Additionally num_trim_bases is set to > 0 then the specified number of bases is always trimmed from every read prior to creating the candidate sequences  

If a sequence can be collapsed to multiple longer sequences, then greedily allocate it to the most supported sequence. The total base qual supporting each base is retained for later matching). 

### Identify anchors and call CDR3 sequences  

For merged sequences, find a V and J anchor simply read out the CDR3 sequence between the V-CYS and J-PHE/J-TRP anchors.  

For each V anchored only consensus sequence we search for candidate J-PHE/J-TRP anchor sequences. To do this we compare each complete 10 amino acid kmer downstream of the V-CYS sequence to the set of known 10 amino acid J anchor sequences by summing the log likelihoods from the BLOSUM62 substitution matrix. Truncated partial anchor sequences of 1 or more amino acids are also checked for the final 1-9 amino acids of the consensus sequence. A similarity score is calculated for each anchor sequence as follows: 

Similarity Score = 3 * Amino Acid length -6 - SUM[Self BLOSUM62 - Anchor BLOSUM62]] 

If the max similarity score to any anchor sequence is greater than 0 we deem it to be a CDR3 sequence. If 2 candidate anchors share the same score, then rank first by inframe and then by whichever CDR3 sequence is closest to 13 amino acids in length. For each J anchored only read/fragment we similarly search for candidate V-CYS anchor sequences in the consensus sequence 

### Collapse consensus sequences 

A sequence is collapsed into another sequence if it is identical at all bases with high quality support across both anchors and candidate CDR3 sequence. Sequences may also be collapsed where they have a single base difference with up to 2 reads of high quality support.  

Each V only anchored read is also checked for partial overlap with each J only anchored read. If they exactly match with more than 15 nucleotides of exact match (allowing for low quality bases) then they are also collapsed into a single sequence. 


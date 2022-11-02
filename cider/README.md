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


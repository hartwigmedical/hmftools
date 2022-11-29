import pandas as pd

import argparse
import os, sys, time
import re
import tempfile
import subprocess

import logging

logger = logging.getLogger(__name__)

logging.basicConfig(stream=sys.stdout,
                    format='%(asctime)s.%(msecs)03d | %(threadName)s | %(levelname)s | %(name)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.DEBUG)

GENE_REGION_TOLERANCE = 50
FLANKING_BASES = 50

'''
   	    qseqid means Query Seq-id
   	       qgi means Query GI
   	      qacc means Query accesion
   	   qaccver means Query accesion.version
   	      qlen means Query sequence length
   	    sseqid means Subject Seq-id
   	 sallseqid means All subject Seq-id(s), separated by a ';'
   	       sgi means Subject GI
   	    sallgi means All subject GIs
   	      sacc means Subject accession
   	   saccver means Subject accession.version
   	   sallacc means All subject accessions
   	      slen means Subject sequence length
   	    qstart means Start of alignment in query
   	      qend means End of alignment in query
   	    sstart means Start of alignment in subject
   	      send means End of alignment in subject
   	      qseq means Aligned part of query sequence
   	      sseq means Aligned part of subject sequence
   	    evalue means Expect value
   	  bitscore means Bit score
   	     score means Raw score
   	    length means Alignment length
   	    pident means Percentage of identical matches
   	    nident means Number of identical matches
   	  mismatch means Number of mismatches
   	  positive means Number of positive-scoring matches
   	   gapopen means Number of gap openings
   	      gaps means Total number of gaps
   	      ppos means Percentage of positive-scoring matches
   	    frames means Query and subject frames separated by a '/'
   	    qframe means Query frame
   	    sframe means Subject frame
   	      btop means Blast traceback operations (BTOP)
   	    staxid means Subject Taxonomy ID
   	  ssciname means Subject Scientific Name
   	  scomname means Subject Common Name
   	sblastname means Subject Blast Name
   	 sskingdom means Subject Super Kingdom
   	   staxids means unique Subject Taxonomy ID(s), separated by a ';'
   			 (in numerical order)
   	 sscinames means unique Subject Scientific Name(s), separated by a ';'
   	 scomnames means unique Subject Common Name(s), separated by a ';'
   	sblastnames means unique Subject Blast Name(s), separated by a ';'
   			 (in alphabetical order)
   	sskingdoms means unique Subject Super Kingdom(s), separated by a ';'
   			 (in alphabetical order) 
   	    stitle means Subject Title
   	salltitles means All Subject Title(s), separated by a '<>'
   	   sstrand means Subject Strand
   	     qcovs means Query Coverage Per Subject
   	   qcovhsp means Query Coverage Per HSP
   	    qcovus means Query Coverage Per Unique Subject (blastn only)
'''

columns = ["qseqid", "qlen", "sseqid", "stitle", "pident", "qcovs", "length", "mismatch", "gapopen",
           "qstart", "qend", "sstart", "send", "qframe", "sframe", "evalue",
           "bitscore", "qseq", "sseq"]

# blastn -db GCF_000001405.39_top_level -query query_seq.fa -task blastn-short -outfmt \
# "6 qseqid sseqid stitle pident qcovs length mismatch gapopen qstart qend sstart send qframe sframe frames evalue bitscore qseq sseq"
def run_blastn(blastn_path, fa_file, blastn_out_file, threads):
    logger.info(f"running blastn({blastn_path}) on {fa_file}")

    cmd = [blastn_path,
           "-db", "GCF_000001405.39_top_level",
           "-query", fa_file,
           "-task", "blastn",
           "-outfmt", "10 " + " ".join(columns),
           "-evalue", "1",
           "-mt_mode", "0",
           "-num_threads", str(threads),
           "-out", blastn_out_file]

    proc = subprocess.run(cmd, cwd=".")
    if proc.returncode == 0:
        logger.info("blastn success")
    else:
        logger.error(f"blastn failed")

class Gene:
    def __init__(self, gene_name, chromosome: str, strand: int, gene_start, gene_end, vdj_type):
        self.gene_name = gene_name
        self.chromosome = str(chromosome)
        self.strand = strand
        self.gene_start = gene_start
        self.gene_end = gene_end
        self.vdj_type = vdj_type

def get_ensembl(ensembl_gene_data, only_vdj):
    ensembl_df = pd.read_csv(ensembl_gene_data)

    # manually add in the KDE section
    # v37: IGKKDE	01	2	89131735	89132285    -
    # v38: IGKKDE	01	chr2	88832222	88832772	-
    ensembl_df = ensembl_df.append(
        {'GeneId': None, 'GeneName': 'KDE', "Chromosome": "2", "Strand": -1, "GeneStart": 89131735,
         "GeneEnd": 89132285}, ignore_index=True)

    ensembl_df = ensembl_df.append(
        {'GeneId': None, 'GeneName': 'KDE', "Chromosome": "2", "Strand": -1, "GeneStart": 88832222,
         "GeneEnd": 88832772}, ignore_index=True)

    # KINT, obtained by BLAT of CACCGCGCTCTTGGGGCAGCCGCCTTGCCGCTAGTGGCCGTGGCCACCCTGTGTCTGCCCGATT
    # v37
    ensembl_df = ensembl_df.append(
        {'GeneId': None, 'GeneName': 'KINT', "Chromosome": "2", "Strand": -1, "GeneStart": 89159145,
         "GeneEnd": 89159209}, ignore_index=True)

    # v38
    ensembl_df = ensembl_df.append(
        {'GeneId': None, 'GeneName': 'KINT', "Chromosome": "2", "Strand": -1, "GeneStart": 88859633,
         "GeneEnd": 88859697}, ignore_index=True)

    genes = []
    for idx, row in ensembl_df.iterrows():
        gene_name = row["GeneName"]
        vdj_type = None
        if (gene_name.startswith("IGH") or
            gene_name.startswith("IGK") or
            gene_name.startswith("IGL") or
            gene_name.startswith("TRA") or
            gene_name.startswith("TRB") or
            gene_name.startswith("TRD") or
            gene_name.startswith("TRG")):

            if len(gene_name) > 3 and (gene_name[3] in "VDJC"):
                vdj_type = gene_name[3]

        if gene_name == "KDE":
            vdj_type = "J"

        if gene_name == "KINT":
            vdj_type = "V"

        if not only_vdj or vdj_type:
            genes.append(Gene(row["GeneName"], row["Chromosome"], row["Strand"], row["GeneStart"], row["GeneEnd"], vdj_type))

    return genes

class GeneFinder:
    def __init__(self, genes: [Gene]):
        self.genes = genes
        self.chromosome_genes = {}
        for g in genes:
            self.chromosome_genes.setdefault((g.chromosome, g.strand), []).append(g)

    def find_gene(self, chromosome: str, strand: int, start, end):
        matching_genes = []
        k = (chromosome, strand)
        if k in self.chromosome_genes:
            gene_list = self.chromosome_genes[k]
            for g in gene_list:
                assert(g.chromosome == chromosome)
                assert(g.strand == strand)
                if g.gene_start - GENE_REGION_TOLERANCE <= end and g.gene_end + GENE_REGION_TOLERANCE >= start:
                    matching_genes.append(g)
        return matching_genes

def process_blastn_result(blastn_csv, gene_finder):

    df = pd.read_csv(blastn_csv,
                     names=["qseqid", "qlen", "sseqid", "stitle_chr", "stitle_assembly", "pident", "qcovs", "length",
                            "mismatch", "gapopen",
                            "qstart", "qend", "sstart", "send", "qframe", "sframe", "evalue",
                            "bitscore", "qseq", "sseq"])

    # add chromosome number if it applies
    def parse_chromosome(stitle):
        m = re.match("^Homo sapiens chromosome (\w+)$", stitle)
        if m:
            return m.group(1)
        return None

    df["chromosome"] = df["stitle_chr"].apply(lambda x: parse_chromosome(x))

    def identify_gene(row):
        sstart = row["sstart"]
        send = row["send"]
        # in negative strand, start and end would be reversed
        if sstart > send:
            sstart, send = send, sstart
        matching_genes = gene_finder.find_gene(row["chromosome"], row["sframe"], sstart, send)

        if not matching_genes:
            return None

        if len(matching_genes) == 1:
            return matching_genes[0]

        # when there are multiple matches, we see if any of them matches IGKJ IGLV etc
        for gene in matching_genes:
            if gene.vdj_type is not None:
                return gene

        return matching_genes[0]

    #start = time.time()
    #print("getting gene info")
    gene_col = df.apply(identify_gene, axis=1)
    #print(f"finished getting gene info, time taken: {round(time.time() - start, 2)}s")

    df["gene"] = gene_col.apply(lambda x: x.gene_name if x else None)
    df["gene_type"] = gene_col.apply(lambda x: x.vdj_type if x else None)
    df["gene_start"] = gene_col.apply(lambda x: x.gene_start if x else None)
    df["gene_end"] = gene_col.apply(lambda x: x.gene_end if x else None)

    aligns_to_keep = df.groupby("qseqid").apply(choose_alignments)
    # convert from series of list into just a flat serie
    idx_to_keep = aligns_to_keep.explode()
    filtered_df = df.filter(items=idx_to_keep, axis=0)
    del df

    def get_gene_type_df(df, gene_type):
        seq_with_gene_type = df[df["gene_type"] == gene_type].groupby("qseqid")[["gene", "pident", "qstart", "qend"]].first()
        seq_with_gene_type.columns = [gene_type.lower() + c.capitalize() for c in seq_with_gene_type.columns]
        return seq_with_gene_type.reset_index()

    # now we have to combine the data and mark
    seq_with_v = get_gene_type_df(filtered_df, "V")
    seq_with_d = get_gene_type_df(filtered_df, "D")
    seq_with_j = get_gene_type_df(filtered_df, "J")

    # merge together
    seq_gene_df = seq_with_v.merge(seq_with_d, how="outer", on="qseqid").merge(seq_with_j, how="outer", on="qseqid")

    # fine genes that matches ref
    seq_match_ref = pd.DataFrame(filtered_df.groupby("qseqid").apply(lambda x: ((x["qlen"] - x["length"]) < 10).any()),
                                 columns=["blastnFullMatch"]).reset_index()
    seq_gene_df = seq_gene_df.merge(seq_match_ref, how="outer", on="qseqid")

    def blastn_status(row):
        if row["blastnFullMatch"]:
            return "NO_REARRANGEMENT"
        if row["vGene"]:
            if row["dGene"]:
                if row["jGene"]:
                    return "V_D_J"
                else:
                    return "V_D"
            elif row["jGene"]:
                return "V_J"
            else:
                return "V_ONLY"
        else:
            if row["dGene"]:
                if row["jGene"]:
                    return "D_J"
                else:
                    return "D_ONLY"
            elif row["jGene"]:
                return "J_ONLY"
            else:
                return "NO_VDJ_ALIGNMENT"

    seq_gene_df["blastnStatus"] = seq_gene_df.fillna("").apply(blastn_status, axis=1).fillna("NO_ALIGNMENT")

    # now we can merge this back to the original sequence
    return seq_gene_df

def write_fasta(vdj_df, fasta_filename, filter_cdr3):
    if filter_cdr3:
        vdj_df = vdj_df[vdj_df["cdr3AA"] == filter_cdr3]
    with open(fasta_filename, "w") as f:
        for idx, row in vdj_df.iterrows():
            if "DUPLICATE" not in row['filter'] and "MATCHES_REF" not in row['filter']:
                full_seq = row["fullSeq"]
                vdj_seq = row["vdjSeq"]
                i = full_seq.index(vdj_seq)
                query_seq = full_seq[max(i - FLANKING_BASES, 0): min(i + len(vdj_seq) + FLANKING_BASES, len(full_seq))]
                f.write(f">{row['cdr3AA']}\n")
                f.write(f"{query_seq}\n")
                #logger.debug(f"i: {i}, vdj_seq: {vdj_seq}, query_seq: {query_seq}, full_seq: {full_seq}")

def main():
    start = time.time()
    logger.info(f"starting cider blastn")
    parser = argparse.ArgumentParser(description='Running blastn on CIDER sequences')
    parser.add_argument('--in_tsv', help='input vdj_seq tsv file', required=True)
    parser.add_argument('--out_tsv', help='output vdj_seq tsv file', required=True)
    parser.add_argument('--temp_dir', help='blastn working directory')
    parser.add_argument('--blastn', help='path to the blastn binary', required=True)
    parser.add_argument('--ensembl', help='path to the ensembl_gene_data csv', required=True)
    parser.add_argument('--threads', type=int, default=1, help='number of threads to use in blastn')
    parser.add_argument('--filter_cdr3', help='only run on this CDR3')
    args = parser.parse_args()

    # read in all the VDJ sequences and blast them
    vdj_df = pd.read_csv(args.in_tsv, sep="\t")

    file_prefix = os.path.basename(args.in_tsv).split('.', maxsplit=1)[0]

    # find all the sequences and put them into a fasta file
    #fd, fa_file = tempfile.mkstemp(suffix=".fa")
    #with os.fdopen(fd, "w") as f:
    fa_file = f"{args.temp_dir}/{file_prefix}.blastn.fa"
    write_fasta(vdj_df, fa_file, args.filter_cdr3)
    blastn_csv = f"{args.temp_dir}/{file_prefix}.blastn.csv"

    # now run blastn on all of them
    run_blastn(args.blastn, fa_file, blastn_csv, args.threads)

    gene_finder = GeneFinder(get_ensembl(args.ensembl, True))
    seq_gene_df = process_blastn_result(blastn_csv, gene_finder)

    # now merge the results back into the vdj_df
    vdj_df = vdj_df.merge(seq_gene_df, how="left", left_on="cdr3AA", right_on="qseqid").drop("qseqid", axis=1)
    vdj_df["blastnStatus"] = vdj_df["blastnStatus"].fillna("SKIPPED_BLASTN")
    vdj_df.to_csv(args.out_tsv, sep="\t")

    elapsed_sec = time.time() - start
    minute = int(elapsed_sec / 60)
    sec = round(elapsed_sec % 60)
    logger.info(f"finished processing {args.in_tsv}, time taken: {minute}m {sec}s")

# we want to choose alignments, based on following:
# 1. if there is one alignment that can encompass the whole range we will take it
#    as it indicates that this sequence matches ref genome
# 2. then for each aligned section we choose the one that is aligned V, D, J or KDE
#    genes. If there are more than one, we choose the highest score one
# 3. next we just choose whatever is left, removing duplicates
def choose_alignments(df):
    rows_to_keep = []

    for idx, row in df.iterrows():
        if row["qlen"] - row["length"] < 5:
            # this row aligns with the whole sequence
            return [idx]

        # check against previous ones and see if can replace any
        i = 0
        while True:
            if i == len(rows_to_keep):
                rows_to_keep.append((idx, row))
                break
            _, r = rows_to_keep[i]
            same_range = abs(r["qstart"] - row["qstart"]) <= 15 and abs(r["qend"] - row["qend"]) <= 15
            contain = ((r["qstart"] <= row["qstart"] and r["qend"] >= row["qend"]) or
                       (r["qstart"] >= row["qstart"] and r["qend"] <= row["qend"]))
            if same_range or contain:
                # these are essentially same ranges
                if row["gene_type"] and not r["gene_type"]:
                    # existing one is not a VDJ gene, replace
                    rows_to_keep[i] = (idx, row)
                elif (row["gene_type"] is None) == (r["gene_type"] is None) and row["length"] > r["length"]:
                    # choose longer one
                    rows_to_keep[i] = (idx, row)
                break
            i += 1

    # next we want to work out which rows are actually duplicates of others
    return [idx for idx, row in rows_to_keep]

def unittest():
    # first test that if there is one range that encampass them all, we take that one
    #         qlen, length, qstart, qend, gene_type
    data = {1: [90,   30,      1,    30,     "V"],
            2: [90,   60,      31,   90,     "J"],
            3: [90,   88,      1,    88,     None]}
    df = pd.DataFrame.from_dict(data, orient='index', columns=["qlen", "length", "qstart", "qend", "gene_type"])
    assert(len(choose_alignments(df)) == 1)
    assert (choose_alignments(df)[0] == 3)

    # when there are multiple ones to choose from, we choose the one with V / J genes
    #         qlen, length, qstart, qend, gene_type
    data = {1: [90,   30,      1,    30,     None],
            2: [90,   27,      3,    29,     "V"],
            3: [90,   28,      3,    30,     "V"],
            4: [90,   59,      31,   90,     "J"],
            5: [90,   58,      31,   90,     "J"],
            6: [90,   60,      31,   90,     None]}
    df = pd.DataFrame.from_dict(data, orient='index', columns=["qlen", "length", "qstart", "qend", "gene_type"])
    assert(len(choose_alignments(df)) == 2)
    # we should keep row 2 and 4
    assert (choose_alignments(df)[0] == 3)
    assert (choose_alignments(df)[1] == 4)

if __name__ == "__main__":
    main()
    #unittest()

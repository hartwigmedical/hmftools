import pandas as pd

import argparse
import os, sys, time
import re
import tempfile
import subprocess
import unittest

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

    logger.info(cmd)

    proc = subprocess.run(cmd, cwd=".")
    if proc.returncode == 0:
        logger.info("blastn success")
    else:
        logger.error(f"blastn failed")

    # run again with normal output
    cmd = [blastn_path,
           "-db", "GCF_000001405.39_top_level",
           "-query", fa_file,
           "-task", "blastn",
           "-evalue", "1",
           "-mt_mode", "0",
           "-num_threads", str(threads),
           "-out", blastn_out_file + ".human_readable"]

    logger.info(cmd)

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

        if not vdj_type:
            continue

        # remove pseudo genes
        if gene_name.endswith("P") or "OR" in gene_name:
            continue

        # also remove Roman numeral ones, like IGHVII, they are pseudo genes
        if len(gene_name) > 4 and gene_name[4].isalpha():
            continue

        if not only_vdj or vdj_type:
            genes.append(Gene(row["GeneName"], row["Chromosome"], row["Strand"], row["GeneStart"], row["GeneEnd"], vdj_type))

    # manually add in the KDE section
    # v38: IGKDEL	01	chr2	88832222	88833331	-
    genes.append(Gene(gene_name='IGKDEL', chromosome="2", strand=-1,
                      gene_start=88832222, gene_end=88833331, vdj_type="J"))
    
    # KINTR, obtained by BLAT of CACCGCGCTCTTGGGGCAGCCGCCTTGCCGCTAGTGGCCGTGGCCACCCTGTGTCTGCCCGATT
    # v38
    genes.append(Gene(gene_name='IGKINTR', chromosome="2", strand=-1,
                      gene_start=88859633, gene_end=88859697, vdj_type="V"))

    #for g in genes:
    #    print(f"Gene(gene_name='{g.gene_name}', chromosome='{g.chromosome}', strand={g.strand}, gene_start={g.gene_start}, gene_end={g.gene_end}, vdj_type='{g.vdj_type}')")

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

def add_gene_data(df, gene_finder):
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
        assert(sstart < send)
        matching_genes = gene_finder.find_gene(row["chromosome"], row["sframe"], sstart, send)

        if not matching_genes:
            return None

        if len(matching_genes) == 1:
            return matching_genes[0]

        # when there are multiple matches, we see if any of them matches IGKJ IGLV etc
        # we return the first match, which would be the one with the highest alignment score
        for gene in matching_genes:
            if gene.vdj_type is not None:
                return gene

        return matching_genes[0]

    #start = time.time()
    #print("getting gene info")

    # use the ensembl data to match the alignments to genes
    gene_col = df.apply(identify_gene, axis=1)
    #print(f"finished getting gene info, time taken: {round(time.time() - start, 2)}s")

    df["gene"] = gene_col.apply(lambda x: x.gene_name if x else None)
    df["gene_type"] = gene_col.apply(lambda x: x.vdj_type if x else None)
    df["gene_chr"] = gene_col.apply(lambda x: x.chromosome if x else None)
    df["gene_start"] = gene_col.apply(lambda x: x.gene_start if x else None)
    df["gene_end"] = gene_col.apply(lambda x: x.gene_end if x else None)

# locus must either be the same, or TRA match with TRD
def compatible_locus(gene_segment1, gene_segment2):
    locus1 = gene_segment1[:3]
    locus2 = gene_segment2[:3]
    if locus1 == locus2:
        return True
    if locus1 == "TRA" and locus2 == "TRD":
        return True
    if locus2 == "TRA" and locus1 == "TRD":
        return True
    return False

# from the blastn dataframe, we get the V, D and J alignments
def get_vdj_alignments(df):

    def get_alignment_df(df, vdj_type):
        # vdj is either V, D or J
        # we take the first one of each type
        gene_type_alignment = df[df["gene_type"] == vdj_type].groupby("qseqid")[["gene", "pident", "qstart", "qend"]].first()
        # rename columns
        gene_type_alignment.rename(columns={"gene": vdj_type.lower() + "Gene",
                                            "pident": vdj_type.lower() + "PIdent",
                                            "qstart": vdj_type.lower() + "AlignStart",
                                            "qend": vdj_type.lower() + "AlignEnd"}, inplace=True)
        return gene_type_alignment.reset_index()

    # now we have to combine the data and mark
    v_alignment = get_alignment_df(df, "V")
    d_alignment = get_alignment_df(df, "D")
    j_alignment = get_alignment_df(df, "J")

    return v_alignment, d_alignment, j_alignment

# for each sequence we assign a locus (IGH, TRA etc), by choosing the highest scoring alignment to
# any VDJ gene segment. The locus is then the locus of that alignment
# We then filter out any alignment that are VDJ segments but not compatible to that locus
def filter_by_locus(df: pd.DataFrame) -> pd.DataFrame:

    # create a df of qseqid, locus
    locus_df = df[df["gene_type"].notna()].sort_values("bitscore", ascending=False).groupby("qseqid")[["gene"]].first().reset_index()
    locus_df["locus"] = locus_df["gene"].str.slice(stop=3)
    locus_df.drop(columns="gene", inplace=True)

    # now put the locus into every row
    df = df.merge(locus_df, on="qseqid", how="left")

    # remove any alignment that do not agree with the chosen locus (like IGKV when locus in TRG)
    # we however keep any alignment that do not have ig/tcr gene, so that we can find if sequence aligns somewhere
    # other than the ig/tcr regions
    df = df[df.apply(
        lambda x: pd.isna(x["gene"]) or pd.isna(x["locus"]) or compatible_locus(x["gene"], x["locus"]),
        axis=1)].reset_index(drop=True)

    return df

def process_blastn_result(df: pd.DataFrame, gene_finder: GeneFinder) -> pd.DataFrame:

    add_gene_data(df, gene_finder)
    df = filter_by_locus(df)

    aligns_to_keep = df.groupby("qseqid").apply(choose_alignments)
    # convert from series of list into just a flat series
    idx_to_keep = aligns_to_keep.explode()
    df = df.filter(items=idx_to_keep, axis=0)

    # get the V, D, J alignments as dataframes
    v_alignment, d_alignment, j_alignment = get_vdj_alignments(df)

    # merge the V, D and J genes together
    seq_gene_df = v_alignment.merge(d_alignment, how="outer", on="qseqid").merge(j_alignment, how="outer", on="qseqid")

    # fine genes that matches ref
    seq_match_ref = pd.DataFrame(df.groupby("qseqid").apply(lambda x: ((x["qlen"] - x["length"]) < 10).any()),
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

    seq_gene_df["alignmentStatus"] = seq_gene_df.fillna("").apply(blastn_status, axis=1).fillna("NO_ALIGNMENT")

    # now we can merge this back to the original sequence
    return seq_gene_df

def write_fasta(vdj_df, fasta_filename, filter_cdr3):

    # to allow sequences to match properly, we need to add 50 flanking bases to the VDJ sequences
    # we add those as columns
    def blastn_query_range(row):
        if "DUPLICATE" not in row['filter'] and "MATCHES_REF" not in row['filter']:
            full_seq = row["fullSeq"]
            vdj_seq = row["vdjSeq"]
            i = full_seq.index(vdj_seq)
            return [max(i - FLANKING_BASES, 0), min(i + len(vdj_seq) + FLANKING_BASES, len(full_seq))]
        return [-1, -1]

    # get the query sequence range
    vdj_df[["blastnQueryStart", "blastnQueryEnd"]] = vdj_df.apply(blastn_query_range, axis=1, result_type="expand")

    if filter_cdr3:
        vdj_df = vdj_df[vdj_df["cdr3AA"] == filter_cdr3]

    with open(fasta_filename, "w") as f:
        for idx, row in vdj_df.iterrows():
            if  row['blastnQueryStart'] != -1 and row['blastnQueryEnd'] != -1:
                # start = row['blastnQueryStart']
                # end = row['blastnQueryEnd']
                # logger.info(f"start: {start}, end: {end}")
                query_seq = row["fullSeq"][row['blastnQueryStart']: row['blastnQueryEnd']]
                # use the index as the key in the fasta so we can identify it later
                f.write(f">{idx}\n")
                f.write(f"{query_seq}\n")
                #logger.debug(f"i: {i}, vdj_seq: {vdj_seq}, query_seq: {query_seq}, full_seq: {full_seq}")

def main():
    start = time.time()
    logger.info(f"starting cider blastn")
    parser = argparse.ArgumentParser(description='Running blastn on CIDER sequences')
    parser.add_argument('--in_tsv', help='input vdj_seq tsv file', required=True)
    parser.add_argument('--out_tsv', help='output vdj_seq tsv file', required=True)
    parser.add_argument('--temp_dir', help='directory to create the working files. Uses system temp if not provided')
    parser.add_argument('--blastn', help='path to the blastn binary', required=True)
    parser.add_argument('--ensembl', help='path to the ensembl_gene_data csv', required=True)
    parser.add_argument('--threads', type=int, default=1, help='number of threads to use in blastn')
    parser.add_argument('--filter_cdr3', help='only run on this CDR3')
    args = parser.parse_args()

    # read in all the VDJ sequences and blast them
    vdj_df = pd.read_csv(args.in_tsv, sep="\t")

    vdj_df.drop(columns=[
        "vGene",
        "vPIdent",
        "vAlignStart",
        "vAlignEnd",
        "dGene",
        "dPIdent",
        "dAlignStart",
        "dAlignEnd",
        "jGene",
        "jPIdent",
        "jAlignStart",
        "jAlignEnd",
        "alignmentStatus"], inplace=True, errors='ignore')

    # strip out the tsv / csv / .gz
    file_prefix = re.sub("(\.[tc]sv\.gz)|(\.[tc]sv)", "", os.path.basename(args.in_tsv))

    if args.temp_dir:
        fa_file = f"{args.temp_dir}/{file_prefix}.blastn.fa"
        blastn_csv = f"{args.temp_dir}/{file_prefix}.blastn.csv"
    else:
        fa_file = tempfile.NamedTemporaryFile(delete=False, suffix=".blastn.fa").name
        blastn_csv = tempfile.NamedTemporaryFile(delete=False, suffix=".blastn.csv").name

    # find all the sequences and put them into a fasta file
    write_fasta(vdj_df, fa_file, args.filter_cdr3)

    # now run blastn on all of them
    run_blastn(args.blastn, fa_file, blastn_csv, args.threads)

    # when reading, the stitle inside contains comma, and we want to split that into chr and assemply
    columns_in = columns.copy()
    stitle_index = columns_in.index("stitle")
    del columns_in[stitle_index]
    columns_in.insert(stitle_index, "stitle_chr")
    columns_in.insert(stitle_index + 1, "stitle_assemply")
    blastn_df = pd.read_csv(blastn_csv, names=columns_in)

    gene_finder = GeneFinder(get_ensembl(args.ensembl, True))
    seq_gene_df = process_blastn_result(blastn_df, gene_finder)

    if not args.temp_dir:
        os.unlink(fa_file)
        os.unlink(blastn_csv)

    # now merge the results back into the vdj_df
    vdj_df = vdj_df.merge(seq_gene_df, how="left", left_index=True, right_on="qseqid").drop("qseqid", axis=1)
    vdj_df["alignmentStatus"] = vdj_df["alignmentStatus"].fillna("SKIPPED_ALIGN")

    # fix the V D J aligned positions to be in terms of the full seq
    for s in ['v', 'd', 'j']:
        vdj_df[f"{s}AlignStart"] = vdj_df[f"{s}AlignStart"] + vdj_df["blastnQueryStart"]
        vdj_df[f"{s}AlignEnd"] = vdj_df[f"{s}AlignEnd"] + vdj_df["blastnQueryEnd"]

    vdj_df.drop(columns=["blastnQueryStart", "blastnQueryEnd"], inplace=True)

    vdj_df.to_csv(args.out_tsv, sep="\t", float_format='%.10g', index=False)

    elapsed_sec = time.time() - start
    minute = int(elapsed_sec / 60)
    sec = round(elapsed_sec % 60)
    logger.info(f"finished processing {args.in_tsv}, time taken: {minute}m {sec}s")

# This function operates on a dataframe that is already grouped by qseqid
# we want to choose alignments, based on following:
# 1. if there is one alignment that can encompass the whole range we will take it
#    as it indicates that this sequence matches ref genome
# 2. then for each aligned section we choose the one that is aligned V, D, J or KDE
#    genes. If there are more than one, we choose the highest score one
# 3. next we just choose whatever is left, removing duplicates
def choose_alignments(df: pd.DataFrame):
    # make sure we already grouped by qseqid
    assert(len(df["qseqid"].unique()) == 1)

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
                elif ((row["gene_type"] is None) == (r["gene_type"] is None) and
                      row["bitscore"] > r["bitscore"]):
                    # choose higher scoring one if they both have same V J type
                    rows_to_keep[i] = (idx, row)
                break
            i += 1

    # next we want to work out which rows are actually duplicates of others
    return [idx for idx, row in rows_to_keep]

class UnitTests(unittest.TestCase):

    def test_choose_alignment(self):
        # first test that if there is one range that encampass them all, we take that one
        #         qseqid, qlen, length, qstart, qend, gene_type, bitscore
        data = {1: [7, 90,   30,      1,    30,     "V",  120],
                2: [7, 90,   60,      31,   90,     "J",  150],
                3: [7, 90,   88,      1,    88,     None, 250]}
        df = pd.DataFrame.from_dict(data, orient='index', columns=["qseqid", "qlen", "length", "qstart", "qend", "gene_type", "bitscore"])
        self.assertEqual(len(choose_alignments(df)), 1)
        self.assertEqual(choose_alignments(df)[0], 3)

        # when there are multiple ones to choose from, we choose the one with V / J genes
        #         qseqid, qlen, length, qstart, qend, gene_type, bitscore
        data = {1: [7, 90,   30,      1,    30,     None, 200],
                2: [7, 90,   27,      3,    29,     "V",  150],
                3: [7, 90,   28,      3,    30,     "V",  250],
                4: [7, 90,   59,      31,   90,     "J",  100],
                5: [7, 90,   58,      31,   90,     "J",  200],
                6: [7, 90,   60,      31,   90,     None, 100]}
        df = pd.DataFrame.from_dict(data, orient='index', columns=["qseqid", "qlen", "length", "qstart", "qend", "gene_type", "bitscore"])
        self.assertEqual(len(choose_alignments(df)), 2)
        # we should keep row 2 and 5 cause of higher bitscore
        self.assertEqual(choose_alignments(df)[0], 3)
        self.assertEqual(choose_alignments(df)[1], 5)

    def test_gene_finder(self):
        genes = [
            Gene(gene_name='IGHV1-18', chromosome="14", strand=-1, gene_start=106184899, gene_end=106185194, vdj_type='V'),
            Gene(gene_name='IGHJ5', chromosome="14", strand=-1, gene_start=105863814, gene_end=105863864, vdj_type='J'),
            Gene(gene_name='IGHD6-25', chromosome="14", strand=-1, gene_start=105881539, gene_end=105881556, vdj_type='D')
        ]
        gene_finder = GeneFinder(genes)

        # test we can correctly annotate the genes
        # first test that if there is one range that encampass them all, we take that one
        #         qseqid, qlen, stitle, length, qstart, qend, sstart, send, sframe, bitscore
        data = {1: [7, 90, "Homo sapiens chromosome 14", 30, 1, 31,   106184899,     106184910, -1, 100],
                2: [7, 90, "Homo sapiens chromosome 14", 60, 31, 91,  105863814,     105863820, -1, 120],
                3: [7, 90, "Homo sapiens chromosome 14", 88, 1, 89,   105881539,     105881541, -1, 130],
                4: [7, 90, "Homo sapiens chromosome 14", 88, 1, 89,   102000000,     102000090, -1, 150]}
        df = pd.DataFrame.from_dict(data, orient='index',
                                    columns=["qseqid", "qlen", "stitle_chr", "length", "qstart", "qend", "sstart", "send", "sframe", "bitscore"])

        # now add gene data
        add_gene_data(df, gene_finder)
        self.assertEqual(df.iloc[0]["gene"], "IGHV1-18")
        self.assertEqual(df.iloc[1]["gene"], "IGHJ5")
        self.assertEqual(df.iloc[2]["gene"], "IGHD6-25")
        self.assertIsNone(df.iloc[3]["gene"])

    def test_filter_by_locus(self):
        genes = [
            Gene(gene_name='IGHV1-18', chromosome="14", strand=-1, gene_start=106184899, gene_end=106185194, vdj_type='V'),
            Gene(gene_name='IGHJ5', chromosome="14", strand=-1, gene_start=105863814, gene_end=105863864, vdj_type='J'),
            Gene(gene_name='IGHD6-25', chromosome="14", strand=-1, gene_start=105881539, gene_end=105881556, vdj_type='D'),
            Gene(gene_name='IGKJ2', chromosome='2', strand=-1, gene_start=88861525, gene_end=88861563, vdj_type='J'),
            Gene(gene_name='IGKV1-22', chromosome='2', strand=-1, gene_start=89170775, gene_end=89171212, vdj_type='V'),
            Gene(gene_name='TRBV3-1', chromosome='7', strand=1, gene_start=142308542, gene_end=142309048, vdj_type='V'),
            Gene(gene_name='TRBD1', chromosome='7', strand=1, gene_start=142786213, gene_end=142786224, vdj_type='D')
        ]
        gene_finder = GeneFinder(genes)
        #         qseqid, qlen, stitle, length, qstart, qend, sstart, send, sframe, bitscore, pident
        data = {1: [7, 150, "Homo sapiens chromosome 14", 90, 1,  90,   106184899,     106184910, -1, 150, 100], # IGHV1
                2: [7, 150, "Homo sapiens chromosome 7", 90, 1, 90, 142308542, 142309048, 1, 100, 100], # TRBV3
                3: [7, 150, "Homo sapiens chromosome 14", 30, 91, 120,  105863814,     105863820, -1, 120, 100], # IGHJ5
                4: [7, 150, "Homo sapiens chromosome 14", 30, 120, 150, 105881539,     105881541, -1, 130, 99.1], # IGHD6
                5: [7, 150, "Homo sapiens chromosome 7",  30, 120, 150, 142786213,    142786224, 1, 100, 99.1], # TRBD1
                6: [7, 150, "Homo sapiens chromosome 14", 88, 1, 89,   102000000,     102000090, -1, 150, 99]}
        blastn_df = pd.DataFrame.from_dict(data, orient='index',
                                    columns=["qseqid", "qlen", "stitle_chr", "length", "qstart", "qend",
                                             "sstart", "send", "sframe", "bitscore", "pident"])

        add_gene_data(blastn_df, gene_finder)
        df = filter_by_locus(blastn_df)

        # should have removed the TRBV3 and TRBD1 rows
        self.assertEqual(len(df), 4)
        self.assertEqual(len(df[df["gene"].str.slice(stop=3) == "TRB"]), 0)
        self.assertEqual(len(df[df["gene"].str.slice(stop=3) == "IGH"]), 3)

    def test_whole(self):
        genes = [
            Gene(gene_name='IGHV1-18', chromosome="14", strand=-1, gene_start=106184899, gene_end=106185194, vdj_type='V'),
            Gene(gene_name='IGHJ5', chromosome="14", strand=-1, gene_start=105863814, gene_end=105863864, vdj_type='J'),
            Gene(gene_name='IGHD6-25', chromosome="14", strand=-1, gene_start=105881539, gene_end=105881556, vdj_type='D'),
            Gene(gene_name='IGKJ2', chromosome='2', strand=-1, gene_start=88861525, gene_end=88861563, vdj_type='J'),
            Gene(gene_name='IGKV1-22', chromosome='2', strand=-1, gene_start=89170775, gene_end=89171212, vdj_type='V')
        ]
        gene_finder = GeneFinder(genes)

        #         qseqid, qlen, stitle, length, qstart, qend, sstart, send, sframe, bitscore, pident
        data = {
                1: [7, 150, "Homo sapiens chromosome 14", 90, 1,  90,   106184899,     106184910, -1, 100, 100], # IGHV1
                2: [7, 150, "Homo sapiens chromosome 14", 30, 91, 120,  105863814,     105863820, -1, 120, 100], # IGHJ5
                3: [7, 150, "Homo sapiens chromosome 14", 30, 120, 150, 105881539,     105881541, -1, 130, 99.1], # IGHD6
                4: [7, 150, "Homo sapiens chromosome 14", 88, 1, 89,   102000000,     102000090, -1, 150, 99],
                }

        # add a higher score one but for a different sequence mapped to IGKV1
        data[5] = [2, 150, "Homo sapiens chromosome 2", 90, 1,  90,   89170800,     89171200, -1, 200, 100]

        blastn_df = pd.DataFrame.from_dict(data, orient='index',
                                    columns=["qseqid", "qlen", "stitle_chr", "length", "qstart", "qend",
                                             "sstart", "send", "sframe", "bitscore", "pident"])
        # now add gene data
        df = process_blastn_result(blastn_df, gene_finder)
        df_qseqid = df[df["qseqid"] == 7]
        self.assertEqual(df_qseqid.iloc[0]["vGene"], "IGHV1-18")
        self.assertEqual(df_qseqid.iloc[0]["jGene"], "IGHJ5")
        self.assertEqual(df_qseqid.iloc[0]["dGene"], "IGHD6-25")

        df_qseqid = df[df["qseqid"] == 2]
        self.assertEqual(df_qseqid.iloc[0]["vGene"], "IGKV1-22")
        self.assertTrue(pd.isna(df_qseqid.iloc[0]["jGene"]))
        self.assertTrue(pd.isna(df_qseqid.iloc[0]["dGene"]))


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "unittest":
        suite = unittest.TestLoader().loadTestsFromTestCase(UnitTests)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()

import pandas as pd

import argparse
import sys

import logging

logger = logging.getLogger(__name__)

logging.basicConfig(stream=sys.stdout,
                    format='%(asctime)s.%(msecs)03d | %(threadName)s | %(levelname)s | %(name)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.DEBUG)

ig_tcr_genes = ['IGLV', 'IGLJ', 'IGLC', 'TRAV', 'TRDV', 'TRDJ', 'TRDC',
       'TRAJ', 'TRAC', 'IGHA', 'IGHE', 'IGHG', 'IGHD', 'IGHM', 'IGHJ',
       'IGHV', 'IGKC', 'IGKJ', 'IGKV', 'TRGC', 'TRGJ', 'TRGV',
       'TRBV', 'TRBD', 'TRBJ', 'TRBC']

# following genes are either not ig/tcr genes, or are pseudogenes
blacklist = ['IGHMBP2', 'IGHGP', 'IGHEP1', 'IGHEP2', 'IGHJ3P', 'IGHJ2P', 'IGHJ1P', 'TRGV5P']

'''
Filter ensembl gene data to get the IG/TCR genes 
'''
class Gene:
    def __init__(self, gene_name, chromosome: str, strand: int, gene_start, gene_end, vdj_type, is_pseudo_gene):
        self.gene_name = gene_name
        self.chromosome = str(chromosome)
        self.strand = strand
        self.gene_start = gene_start
        self.gene_end = gene_end
        self.vdj_type = vdj_type
        self.is_pseudo_gene = is_pseudo_gene

def get_ensembl(ensembl_gene_data):
    ensembl_df = pd.read_csv(ensembl_gene_data)

    genes = []
    for idx, row in ensembl_df.iterrows():

        assert(row["Strand"] == 1 or row["Strand"] == -1)
        gene_name = row["GeneName"]

        if gene_name in blacklist:
            continue

        if (len(gene_name) < 4 or gene_name[:4] not in ig_tcr_genes):
            continue

        vdj_type = gene_name[3]

        if (gene_name == "IGHG2"):
            logger.info(f"gene {vdj_type}")

        if not vdj_type:
            continue

        is_peudogene = False

        # remove pseudo genes and Orphan genes
        if "OR" in gene_name:
            logger.info(f"removing orphan gene: {gene_name}")
            is_peudogene = True
            continue

        # also remove Roman numeral ones, like IGHVII, they are pseudo genes
        if len(gene_name) > 4 and gene_name[4] in "VIX":
            logger.info(f"removing psuedo gene: {gene_name}")
            is_peudogene = True
            continue

        genes.append(
            Gene(row["GeneName"], row["Chromosome"], row["Strand"], row["GeneStart"], row["GeneEnd"], vdj_type, is_peudogene))

    # manually add in the KDE section
    # v38: IGKDEL	01	chr2	88832222	88833331	-
    genes.append(Gene(gene_name='IGKDEL', chromosome="2", strand=-1,
                      gene_start=88832222, gene_end=88833331, vdj_type="J", is_pseudo_gene=False))

    # KINTR, obtained by BLAT of CACCGCGCTCTTGGGGCAGCCGCCTTGCCGCTAGTGGCCGTGGCCACCCTGTGTCTGCCCGATT
    # v38
    genes.append(Gene(gene_name='IGKINTR', chromosome="2", strand=-1,
                      gene_start=88859633, gene_end=88859697, vdj_type="V", is_pseudo_gene=False))

    # for g in genes:
    #    print(f"Gene(gene_name='{g.gene_name}', chromosome='{g.chromosome}', strand={g.strand}, gene_start={g.gene_start}, gene_end={g.gene_end}, vdj_type='{g.vdj_type}')")

    return genes

def write_gene_data(out_tsv: str, genes: [Gene]):
    # make it into object
    cols = {
        'gene': [g.gene_name for g in genes],
        'chromosome': [g.chromosome for g in genes],
        'strand': ['+' if g.strand == 1 else '-' for g in genes],
        'start': [g.gene_start for g in genes],
        'end': [g.gene_end for g in genes],
        'geneSegmentType': [g.vdj_type for g in genes]
        #'psuedoGene': [g.is_pseudo_gene for g in genes]
    }

    df = pd.DataFrame(data=cols)
    df.to_csv(out_tsv, sep="\t", index=False)
    return

def main():
    logger.info(f"starting")
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_tsv', help='output gene tsv file', required=True)
    parser.add_argument('--ensembl', help='path to the ensembl_gene_data csv', required=True)
    args = parser.parse_args()

    genes = get_ensembl(args.ensembl)
    write_gene_data(args.out_tsv, genes)
    logger.info(f"finished writing {len(genes)} genes")

if __name__ == "__main__":
    main()

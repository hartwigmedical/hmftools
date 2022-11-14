import pandas as pd
import json
import argparse

def isX(chromosome: str):
    return chromosome == "X" or chromosome == "chrX"

def isY(chromosome: str):
    return chromosome == "Y" or chromosome == "chrY"

def load_exon_region_df(csv_path, region_size):
    # for each position we also want to add the name of the exon it belongs to, if any
    exons_df = pd.read_csv(csv_path, names=['chromosome', 'start', 'end', 'exon'], header=None, sep='\t',
                           dtype={"chromosome": str})

    # in order to merge with the regions, we need to cut open the regions to smaller ones
    def exon_to_regions(exon_row):
        # we round the start position to 1k
        start_region = int(exon_row['start'] / region_size) * region_size + 1
        end_region = int(exon_row['end'] / region_size + 1) * region_size + 1
        #print(exon_row.start, exon_row.end)
        return [r for r in range(start_region, end_region, region_size)]

    exon_region_df = exons_df.reset_index(drop=True)
    exon_region_df['position'] = exon_region_df.apply(lambda x: exon_to_regions(x), axis=1)
    exon_region_df = exon_region_df[['chromosome', 'exon', 'position']].explode('position').reset_index(drop=True)

    # remove the duplicate
    exon_region_df.drop_duplicates(subset=['chromosome', 'position'], inplace=True)
    return exon_region_df

def load_cobalt_ratio_df(sample_cfg, assume_diploid):
    with open(sample_cfg, 'r') as f:
        samples = json.load(f)

    cobalt_ratio_dfs = []

    try:
        for sample in samples:
            sample_id = sample["sample_id"]
            targeted_cobalt_ratios = sample["targeted_cobalt_ratios"]

            cols_to_drop = ['referenceReadCount', 'referenceGCRatio', 'referenceGCDiploidRatio']

            panel_df = pd.read_csv(targeted_cobalt_ratios, sep="\t", dtype={"chromosome": str}).drop(columns=cols_to_drop)
            panel_df = panel_df.loc[panel_df['tumorGCRatio'] >= 0]

            if assume_diploid:
                gender = sample["gender"]
                # here we do not have the WGS data to check against. We have to assume the tumorGCRatio is
                # mostly 1.0
                # except in the sex chromosome, where we assume it to be 0.5 for X/Y for male, and nan for
                # Y in female.
                wgs_cobalt_ratios = "N/A"
                wgs_df = panel_df.copy(deep=True)
                wgs_df["gender"] = gender

                if gender == 'MALE':
                    # for male, set all X/Y chromosome tumorGCRatio to 0.5
                    wgs_df['tumorGCRatio'] = wgs_df.apply(
                        lambda x: 0.5 if isX(x["chromosome"]) or isY(x["chromosome"]) else 1.0, axis=1)
                elif gender == 'FEMALE':
                    # for female, set Y chromosome tumorGCRatio to nan
                    wgs_df['tumorGCRatio'] = wgs_df.apply(
                        lambda x: float('nan') if isY(x["chromosome"]) else 1.0, axis=1)
                else:
                    raise RuntimeError(f"unknown gender: {gender}")
            else:
                wgs_cobalt_ratios = sample["wgs_cobalt_ratios"]
                wgs_df = pd.read_csv(wgs_cobalt_ratios, sep="\t", dtype={"chromosome": str}).drop(columns=cols_to_drop)
                wgs_df = wgs_df.loc[wgs_df['tumorGCRatio'] >= 0]

            # merge them together
            combine_df = pd.merge(panel_df, wgs_df, how='left', on=['chromosome', 'position'], suffixes=['_panel', '_wgs'])

            combine_df.insert(0, 'sample_id', sample_id)
            cobalt_ratio_dfs.append(combine_df)

            print(f"loaded sample {sample_id}, targeted cobalt={targeted_cobalt_ratios}, wgs cobalt={wgs_cobalt_ratios}")

    except KeyError as e:
        print(f"error loading {sample_cfg}, missing key: {e}")
        raise e

    # merge all into one
    cobalt_ratio_df = pd.concat(cobalt_ratio_dfs).reset_index(drop=True)
    return cobalt_ratio_df

def chromosome_rank(chromosome):
    if chromosome is str:
        chromosome = chromosome.replace('chr', '').upper()
    if chromosome == 'X':
        return 23
    if chromosome == 'Y':
        return 24
    if chromosome == 'MT':
        return 25
    try:
        return int(chromosome)
    except:
        return 1000

epilog = '''Argument for sample_cfg is a json file that contains information for each sample.
For each sample we need to provide the cobalt ratio outputs generated from both the
WGS bam file and targeted bam file.

Example file input for sample_cfg:
[
    {
        "sample_id": "FR16648805",
        "wgs_cobalt_ratios": "WIDE01010081T.cobalt.ratio.tsv.gz",
        "targeted_cobalt_ratios": "FR16648805.cobalt.ratio.tsv.gz"
    },
    {
        "sample_id": "FR16648808",
        "wgs_cobalt_ratios": "WIDE01010241T.cobalt.ratio.tsv.gz",
        "targeted_cobalt_ratios": "FR16648808.cobalt.ratio.tsv.gz"
    }
]
Example file input for target_region:
chromosome\tstart\tend\texon
chr1\t2556664 2556733\t0_TNFRSF14_CODING
chr1\t2557725 2557834\t1_TNFRSF14_CODING
chr1\t2558342 2558468\t2_TNFRSF14_CODING
chr1\t2559822 2559978\t3_TNFRSF14_CODING'''

def main():

    parser = argparse.ArgumentParser(description='''Calculate normalisation for target enrichment process

This program looks at the cobalt ratios of WGS bams and targeted(panel) bams of a set of
samples, and calculate the relative enrichment of each targeted region.
The result is written into a TSV file that can be used by cobalt to correct for the ratios
of bam files generated using the same target enrichment process.''',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=epilog)
    parser.add_argument('--output', help='Output TSV file', required=True)
    parser.add_argument('--sample_cfg', help='Path to sample config file', required=True)
    parser.add_argument('--target_region', help='CSV file with targeted regions', required=True)
    parser.add_argument('--cobalt_window_size', type=int, default=1000,
                        help='Cobalt genome region window size [default=1000]')
    parser.add_argument('--assume_diploid', type=bool, default=False,
                        help='Assume the genome is diploid instead of using the WGS ratios')
    parser.add_argument('--min_enrichment_ratio', type=float, default=0.1,
                        help='Minimum enrichment ratio, window with enrichment ratios < this will be masked')
    args = parser.parse_args()

    if args.assume_diploid:
        print(f"--assume_diploid=True set. Assuming genome is diploid")

    exon_region_df = load_exon_region_df(args.target_region, args.cobalt_window_size)
    cobalt_ratio_df = load_cobalt_ratio_df(args.sample_cfg, args.assume_diploid)

    # we want to merge in the exons, so we know which region is within the panel
    cobalt_ratio_df = cobalt_ratio_df.merge(exon_region_df[['chromosome', 'position', 'exon']], how='left', on=('chromosome', 'position'))
    cobalt_ratio_df['target_region'] = cobalt_ratio_df['exon'].notna()

    # remove all regions that are not in the exons. Maybe think about including more?
    cobalt_ratio_df = cobalt_ratio_df.loc[cobalt_ratio_df['target_region']].reset_index(drop=True)
    cobalt_ratio_df

    # now want to normalise this by the median of GC ratios of the panel regions
    normalisation_df = cobalt_ratio_df[['sample_id', 'tumorGCRatio_panel']].groupby(['sample_id']).agg(['median'])
    normalisation_df.columns = normalisation_df.columns.to_flat_index()
    normalisation_df.columns = [ c[0] + '_' + c[1] for c in normalisation_df.columns]
    normalisation_df

    # use this normalisation to calculate relative enrichment per region
    cobalt_ratio_df = pd.merge(cobalt_ratio_df, normalisation_df, how='left', on='sample_id')
    cobalt_ratio_df['relativeEnrichment'] = cobalt_ratio_df['tumorGCRatio_panel'] / cobalt_ratio_df['tumorGCRatio_wgs'] / cobalt_ratio_df['tumorGCRatio_panel_median']
    cobalt_ratio_df

    print(f"calculating relative enrichment")

    # calculate enrichment
    enrichment_df = cobalt_ratio_df[['chromosome', 'position', 'relativeEnrichment']].groupby(['chromosome','position']).agg(['median', 'mean', 'std', 'min', 'max'])
    enrichment_df.reset_index(inplace=True)
    # flatten the columns
    enrichment_df.columns = [ c[1] if c[1] else c[0] for c in enrichment_df.columns.to_flat_index() ]

    # now save all these regions into a bed file for cobalt to use
    enrichment_file = enrichment_df[['chromosome', 'position', 'median']].reset_index(drop=True)
    enrichment_file.columns = ['chromosome', 'position', 'relativeEnrichment']

    # set a minimum relative enrichment, enrichment ratio below this will be removed
    enrichment_file['relativeEnrichment'] = enrichment_file['relativeEnrichment'].apply(lambda x: x if x >= args.min_enrichment_ratio else float('nan'))

    # merge in the exons again. We must do this cause some regions might be missing from the cobalt data
    enrichment_file = enrichment_file.merge(exon_region_df, on=['chromosome', 'position'], how='outer')
    enrichment_file.drop(columns=['exon'], inplace=True)
    enrichment_file['position'] = enrichment_file['position'].astype(int)
    enrichment_file['relativeEnrichment'] = enrichment_file['relativeEnrichment'].apply(lambda x: 'NaN' if pd.isna(x) else '{:.4f}'.format(x))
    enrichment_file.sort_values(by='position', inplace=True)
    enrichment_file.sort_values(by='chromosome', key=lambda col: col.map(chromosome_rank), kind='stable', inplace=True)

    print(f"writing {enrichment_file.shape[0]} enrichments to {args.output}")
    enrichment_file.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()

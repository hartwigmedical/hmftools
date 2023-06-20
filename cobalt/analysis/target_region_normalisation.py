import pandas as pd
import json
import argparse
import numpy

import os, sys
script_dir = os.path.dirname(__file__)
if script_dir not in sys.path:
    sys.path.append(script_dir)

import cobalt_calc

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

        start_pos = exon_row['start'] + 1 # standard BED file adjustment
        start_region = int(numpy.floor(start_pos / region_size) * region_size) + 1

        end_pos = exon_row['end']

        if (end_pos % 1000) == 0:
            end_region = end_pos
        else:
            end_region = int(numpy.ceil(end_pos / region_size) * region_size)

        #exon_start = exon_row['start'] + 1 # BED file standard
        #start_region = int(exon_start / region_size) * region_size + 1
        #end_region = int(exon_row['end'] / region_size + 1) * region_size + 1
        #print(exon_row.start, exon_row.end)
        return [r for r in range(start_region, end_region, region_size)]

    exon_region_df = exons_df.reset_index(drop=True)
    exon_region_df['position'] = exon_region_df.apply(lambda x: exon_to_regions(x), axis=1)
    exon_region_df = exon_region_df[['chromosome', 'exon', 'position']].explode('position').reset_index(drop=True)
    exon_region_df['position'] = exon_region_df['position'].astype(int)

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

            panel_df = pd.read_csv(targeted_cobalt_ratios, sep="\t", dtype={"chromosome": str})
            panel_df.drop(columns=['referenceReadCount', 'referenceGCRatio', 'referenceGCDiploidRatio'],
                          inplace=True)

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
                wgs_df = pd.read_csv(wgs_cobalt_ratios, sep="\t", dtype={"chromosome": str})
                wgs_df = wgs_df.loc[wgs_df['tumorGCRatio'] >= 0]

            # merge them together, we only need tumorGCRatio from WGS
            combine_df = pd.merge(panel_df, wgs_df[['chromosome', 'position', 'tumorGCRatio']],
                                  how='left', on=['chromosome', 'position'], suffixes=('_panel', '_wgs'))

            combine_df.insert(0, 'sample_id', sample_id)
            cobalt_ratio_dfs.append(combine_df)

            print(f"loaded sample {sample_id}, targeted cobalt={targeted_cobalt_ratios}, wgs cobalt={wgs_cobalt_ratios}")

    except KeyError as e:
        print(f"error loading {sample_cfg}, missing key: {e}")
        raise e

    # merge all into one
    cobalt_ratio_df = pd.concat(cobalt_ratio_dfs).reset_index(drop=True)

    return cobalt_ratio_df

# add two boolean columns: onTarget and offTarget
def label_on_off_target(cobalt_ratio_df, exon_region_df) -> pd.DataFrame:
    # we want to merge in the exons, so we know which region is within the panel
    cobalt_ratio_df = cobalt_ratio_df.merge(exon_region_df[['chromosome', 'position', 'exon']], how='inner',
                                            on=('chromosome', 'position'))
    # on target ratio
    cobalt_ratio_df["onTarget"] = cobalt_ratio_df['exon'].notna()

    # add another column whether it is included in offtarget ratio
    exclude_from_offtarget_df = exon_region_df[["chromosome", "position"]].reset_index(drop=True)
    exclude_from_offtarget_df["position"] = exclude_from_offtarget_df["position"].apply(
        lambda x: [v for v in range(x - cobalt_calc.OFF_TARGET_MIN_DISANCE_FROM_ON_TARGET,
                                    x + cobalt_calc.OFF_TARGET_MIN_DISANCE_FROM_ON_TARGET + cobalt_calc.WINDOW_SIZE,
                                    cobalt_calc.WINDOW_SIZE)])
    exclude_from_offtarget_df = exclude_from_offtarget_df.explode("position").drop_duplicates()
    exclude_from_offtarget_df["excludeFromOffTarget"] = True

    cobalt_ratio_df = cobalt_ratio_df.merge(exclude_from_offtarget_df, on=["chromosome", "position"], how="left")
    cobalt_ratio_df["offTarget"] = ~(cobalt_ratio_df["excludeFromOffTarget"].fillna(False))
    cobalt_ratio_df.drop(columns="excludeFromOffTarget", inplace=True)

    return cobalt_ratio_df

def chromosome_rank(chromosome):
    if isinstance(chromosome, str):
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

def create_normalisation_df(cobalt_ratio_df, row_mask) -> pd.DataFrame:

    # make a copy with row mask applied
    cobalt_ratio_df = cobalt_ratio_df[row_mask].reset_index(drop=True)

    # calculate a gc normalised value of all targeted regions
    cobalt_ratio_df["tumorGCRatio_panel"] = cobalt_ratio_df.groupby("sample_id", as_index=False) \
        .apply(lambda x: cobalt_calc.gc_normalise(x, "tumorReadCount", cobalt_calc.included_in_median_calc)) \
        .droplevel(0)

    # use this normalisation to calculate relative enrichment per region
    cobalt_ratio_df['relativeEnrichment'] = cobalt_ratio_df['tumorGCRatio_panel'] / cobalt_ratio_df['tumorGCRatio_wgs']

    # cobalt_ratio_file = 'cobalt_ratio_df.tsv'
    # cobalt_ratio_df.to_csv(cobalt_ratio_file, sep="\t", index=False)

    print(f"calculating relative enrichment")

    # calculate enrichment
    normalisation_df = cobalt_ratio_df[['chromosome', 'position', 'relativeEnrichment']] \
        .groupby(['chromosome', 'position']).agg(['median']).reset_index()

    # flatten the columns
    normalisation_df.columns = ['chromosome', 'position', 'relativeEnrichment']
    return normalisation_df

def create_onoff_target_df(cobalt_ratio_df, exon_region_df, gc_profile_df, min_enrichment_ratio) -> pd.DataFrame:
    # merge in gc profiles
    cobalt_ratio_df = cobalt_ratio_df.merge(gc_profile_df, on=["chromosome", "position"], how="left")

    cobalt_ratio_df = label_on_off_target(cobalt_ratio_df, exon_region_df)

    on_target_normalisation_df = create_normalisation_df(cobalt_ratio_df, cobalt_ratio_df["onTarget"])
    #on_target_normalisation_df["onTarget"] = True
    #on_target_normalisation_df["offTarget"] = False
    #off_target_normalisation_df = create_normalisation_df(cobalt_ratio_df, cobalt_ratio_df["offTarget"])
    #off_target_normalisation_df["onTarget"] = False
    #off_target_normalisation_df["offTarget"] = True

    # merge into one
    enrichment_df = on_target_normalisation_df

    # set a minimum relative enrichment, enrichment ratio below this will be removed
    enrichment_df['relativeEnrichment'] = enrichment_df['relativeEnrichment'].apply(lambda x: x if x >= min_enrichment_ratio else float('nan'))
    return enrichment_df

epilog = '''Argument for sample_cfg is a json file that contains information for each sample.
For each sample we need to provide the cobalt ratio outputs generated from both the
WGS bam file and targeted bam file.

Example file input for sample_cfg:
[
    {
        "sample_id": "FR16648805",
        "wgs_cobalt_ratios": "WIDE01010081T.cobalt.ratio.tsv.gz",
        "targeted_cobalt_ratios": "FR16648805.cobalt.ratio.tsv.gz",
        "gender": "MALE"
    },
    {
        "sample_id": "FR16648808",
        "wgs_cobalt_ratios": "WIDE01010241T.cobalt.ratio.tsv.gz",
        "targeted_cobalt_ratios": "FR16648808.cobalt.ratio.tsv.gz",
        "gender": "FEMALE"
    }
]
Argument for target_region is a bed file what describe all genome regions that are captured:
chromosome\tstart\tend\texon
chr1\t2556664 2556733\t0_TNFRSF14_CODING
chr1\t2557725 2557834\t1_TNFRSF14_CODING
chr1\t2558342 2558468\t2_TNFRSF14_CODING
'''

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
    parser.add_argument('--gc_profile', help='TSV file with gc profiles', required=True)
    parser.add_argument('--cobalt_window_size', type=int, default=1000,
                        help='Cobalt genome region window size [default=1000]')
    parser.add_argument('--assume_diploid', action='store_true',
                        help='Assume the genome is diploid instead of using the WGS ratios')
    parser.add_argument('--min_enrichment_ratio', type=float, default=0.1,
                        help='Minimum enrichment ratio, window with enrichment ratios < this will be masked')
    args = parser.parse_args()

    if args.assume_diploid:
        print(f"--assume_diploid=True set. Assuming genome is diploid")

    exon_region_df = load_exon_region_df(args.target_region, args.cobalt_window_size)
    gc_profile_df = cobalt_calc.load_gc_profile_df(args.gc_profile)

    cobalt_ratio_df = load_cobalt_ratio_df(args.sample_cfg, args.assume_diploid)

    enrichment_file = create_onoff_target_df(cobalt_ratio_df, exon_region_df, gc_profile_df, args.min_enrichment_ratio)

    # merge in the exons again. We must do this cause some regions might be missing from the cobalt data
    enrichment_file = enrichment_file.merge(exon_region_df, on=['chromosome', 'position'], how='outer').drop(columns='exon')
    enrichment_file['relativeEnrichment'] = enrichment_file['relativeEnrichment'].apply(lambda x: 'NaN' if pd.isna(x) else '{:.4f}'.format(x))
    enrichment_file.sort_values(by='position', inplace=True)
    enrichment_file.sort_values(by='chromosome', key=lambda col: col.map(chromosome_rank), kind='stable', inplace=True)

    print(f"writing {enrichment_file.shape[0]} enrichments to {args.output}")
    enrichment_file.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()

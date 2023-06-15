import pandas as pd

import argparse
import os, sys, time
from math import floor
from collections import namedtuple

import logging

logger = logging.getLogger(__name__)

logging.basicConfig(stream=sys.stdout,
                    format='%(asctime)s.%(msecs)03d | %(threadName)s | %(levelname)s | %(name)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)

WINDOW_SIZE = 1000
OFF_TARGET_MIN_DISANCE_FROM_ON_TARGET = 2000
MAX_SPARSE_CONSOLIDATE_DISTANCE = 3000000

def load_gc_profile_df(gc_profile_path) -> pd.DataFrame:
    gc_profile_df = pd.read_csv(gc_profile_path,
                             names = ["chromosome", "position", "gcContent", "nonNPercent", "mappablePercent"],
                             sep="\t", dtype={"chromosome": str})

    # gc profile positions are 0, like 2556000, we need to +1 to
    gc_profile_df["position"] = gc_profile_df["position"] + 1

    # drop all the -1 rows
    gc_profile_df = gc_profile_df[gc_profile_df["gcContent"] != -1].reset_index(drop=True)

    # add a gcBucket column.
    # NOTE: cannot use python or panda round function, as it uses a strange banker rounding. We
    # must round 0.5 up to match with Java code
    gc_profile_df["gcBucket"] = ((gc_profile_df["gcContent"] * 100.0) + 0.5).astype(int)
    
    return gc_profile_df


def included_in_median_calc(x):
    if x["tumorReadCount"] < 0:
        return False
    if x["gcBucket"] < 20.0:
        return False
    if x["gcBucket"] > 60.0:
        return False
    if x["mappablePercent"] < (0.85 - 1e-10):
        return False
    if x["chromosome"].endswith("X"):
        return False
    if x["chromosome"].endswith("Y"):
        return False
    return True


def interpolated_median(input: pd.Series) -> float:
    EPSILON = 1e-10

    if input.empty:
        return 0

    values = sorted(input)
    count = len(values)

    median = input.median()

    count_below_med = 0
    count_above_med = 0

    cumulative_count = 0

    for val in values:
        if cumulative_count == count / 2.0 and val != median:
            return median

        if val < (median - EPSILON):
            count_below_med += 1
        elif val > (median + EPSILON):
            count_above_med += 1

        cumulative_count += 1

    l = count_below_med / count
    r = count_above_med / count

    return median - 0.5 + (0.5 - l) / (1 - l - r)


# apply gcratio
def calc_gc_ratio(x: pd.Series, sample_median_read_count, sample_mean_read_count, input_column: str):
    median_normalisation = sample_median_read_count / sample_mean_read_count
    
    if x["gcMedianCount"] == 0 or pd.isna(x["gcMedianCount"]):
        return None
    
    ratio = median_normalisation * x[input_column] / x["gcMedianCount"]
    return ratio

# gc normalise the input column, keep the original df the same
def gc_normalise(df: pd.DataFrame, input_column: str, included_in_median_calc_func):
    # calculate sample median
    # it is non X, Y and mappable regions
    df_median_calc = df[df.apply(included_in_median_calc_func, axis=1)]

    gc_median_df = df_median_calc[["gcBucket", input_column]].groupby("gcBucket").agg([interpolated_median, "count"]).reset_index()
    gc_median_df.columns = ["gcBucket", "gcMedianCount", "windowCount"]

    sample_median_read_count = interpolated_median(df_median_calc[input_column])
    sample_mean_read_count = df_median_calc[input_column].mean()

    #print(f"sample_median_read_count: {sample_median_read_count}, sample_mean_read_count: {sample_mean_read_count}")

    # merge in the gcMedianCount
    # here the index is reset to make sure we get the same index as original df after merge
    df = df.merge(gc_median_df[["gcBucket", "gcMedianCount"]], on="gcBucket", how="left").set_index(df.index)

    return df.apply(lambda x: calc_gc_ratio(x, sample_median_read_count, sample_mean_read_count, input_column), axis=1)

def add_on_target_ratio(df):
    #target_region_gcratio_median = df[df["relativeEnrichment"].notna()]["tumorReadCount"].median()
    #target_region_gcratio_median
    df["onTargetGcRatio"] = df["tumorReadCount"] / df["relativeEnrichment"]
    df["onTargetGcRatio"] = gc_normalise(df, "onTargetGcRatio", lambda x: pd.notna(x["relativeEnrichment"]) and included_in_median_calc(x))

def make_cobalt_dataframe(cobalt_ratio_path, gc_profile_path, targeted_regions_normalisation_path) -> pd.DataFrame:

    gc_profile_df = load_gc_profile_df(gc_profile_path)
    
    # for each position we also want to add the name of the exon it belongs to, if any
    target_regions_normalisation_df = pd.read_csv(targeted_regions_normalisation_path, sep='\t').dropna()

    df = pd.read_csv(cobalt_ratio_path, sep="\t")

    # merge in the gc profile
    df = df.merge(gc_profile_df, on=["chromosome", "position"], how="left")

    # merge in the relative enrichment
    df = df.merge(target_regions_normalisation_df, on=["chromosome", "position"], how="left")

    df["rawTumorGcRatio"] = gc_normalise(df, "tumorReadCount", included_in_median_calc)
    
    add_on_target_ratio(df)

    # add another column whether it is included in offtarget ratio
    exclude_from_offtarget_df = target_regions_normalisation_df[["chromosome", "position"]].reset_index(drop=True)
    exclude_from_offtarget_df["position"] = exclude_from_offtarget_df["position"].apply(
        lambda x: [v for v in range(x - OFF_TARGET_MIN_DISANCE_FROM_ON_TARGET, x + OFF_TARGET_MIN_DISANCE_FROM_ON_TARGET + WINDOW_SIZE, WINDOW_SIZE)])
    exclude_from_offtarget_df = exclude_from_offtarget_df.explode("position").drop_duplicates()
    exclude_from_offtarget_df["excludeFromOffTarget"] = True

    df = df.merge(exclude_from_offtarget_df, on=["chromosome", "position"], how="left")
    df["includedInOffTarget"] = ~(df["excludeFromOffTarget"].fillna(False))
    df.drop(columns="excludeFromOffTarget", inplace=True)

    return df

# drop columns if we want to use this for proper use
def drop_analysis_columns(df):
    df.drop(columns=["includedInOffTarget"], inplace=True)
    df.rename({"onTargetGcRatio": "tumorGCRatio"}, inplace=True)

LowCovBucket = namedtuple('LowCovBucket', ['start', 'end', 'position'])

# assign bucket IDs
# given a list of window positions, return a list of bucket IDs
def assign_bucket_ids(df: pd.DataFrame, consolidation_count: int):

    # window position has to be sorted already, it won't work otherwise
    #assert(window_positions.sort_values())

    # Helper function to round down to window boundary
    def round_down_to_window_boundary(p):
        return floor(p / WINDOW_SIZE) * WINDOW_SIZE + 1

    bucket_ids = [float('nan')] * len(df)

    window_positions = df["position"]
    window_count = 0
    current_bucket_id = 0
    bucket_start = -1
    chromosome = ""

    for i in range(len(df)):
        row = df.iloc[i]
        position = window_positions[i]

        if (chromosome != row["chromosome"] or
            position - bucket_start) >= MAX_SPARSE_CONSOLIDATE_DISTANCE:
            # Do not let the bucket consolidate over 3M bases
            # move to next bucket id
            current_bucket_id += 1

            # Reset bucket start position
            # Also reset window count
            bucket_start = position
            window_count = 0

        if window_count == consolidation_count:
            current_bucket_id += 1

            # Put the bucket boundary in the middle of the two windows
            last_position = window_positions[i - 1]
            bucket_end = round_down_to_window_boundary((last_position + position) * 0.5)

            # Next bucket starts right after this
            # Also reset window count
            bucket_start = bucket_end + WINDOW_SIZE
            window_count = 0

        window_count += 1
        bucket_ids[i] = current_bucket_id

    df["bucket_id"] = bucket_ids

def create_lowcov_bucket(df: pd.DataFrame, consolidate_count):
    # for each chromosome we find the buckets and assign the IDs back
    # df["bucket_id"] = assign_bucket_ids(x, consolidate_count)

    # create a table of the bucket as well
    #df.groupby("bucket_id")["position"].agg(["min", "max"])
    pass


def main():
    #logger.info(f"running cobalt logic")
    parser = argparse.ArgumentParser(description='cobalt logic')
    parser.add_argument('--cobalt_ratio', help='input cobalt ratios tsv file', required=True)
    parser.add_argument('--output', help='output cobalt ratios tsv file', required=True)
    parser.add_argument('--gc_profiles', help='path to the gc profiles file', required=True)
    parser.add_argument('--target_regions_normalisation', help='path to the targeted regions normalisation file', required=True)
    args = parser.parse_args()

    df = make_cobalt_dataframe(args.cobalt_ratio, args.gc_profiles, args.target_regions_normalisation)
    drop_analysis_columns(df)
    df.to_csv(args.output, sep="\t")

if __name__ == "__main__":
    main()

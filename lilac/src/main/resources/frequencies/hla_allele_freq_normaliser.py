import argparse
import logging
logging.basicConfig(level=logging.DEBUG, format="%(asctime)s [%(levelname)-5s] %(message)s", datefmt="%H:%M:%S")

import numpy as np
import pandas as pd
pd.set_option('display.max_columns', 100)
pd.set_option('display.width', 300)


class HlaAlleleFreqNormaliser:

    SAMPLE_SIZE_CAP = 1000
    NOT_PRESENT_FREQUENCY = 0
    MIN_ALLELE_OBSERVATIONS = 2

    def __init__(self, afnd_path: str, output_dir: str):
        self.afnd_path = afnd_path
        self.output_dir = output_dir

        self._frequencies: pd.DataFrame = None
        self._stats: pd.DataFrame = None

    def load_afnd(self) -> None:

        logging.info(f"Loading scraped Allele Frequency Net Database entries from file: {self.afnd_path}")

        df = pd.read_table(self.afnd_path)
        logging.debug(f"Loaded {len(df)} entries")

        logging.info("Parsing column values...")

        df["sample_size"] = df["sample_size"].str.replace(",", "").astype(int)
        df["allele_frequency"] = df["allele_frequency"].astype(str).str.replace("(*)", "").astype(float)
        df["study_weight"] = np.minimum(df["sample_size"], HlaAlleleFreqNormaliser.SAMPLE_SIZE_CAP)
        df["approx_observations"] = df["allele_frequency"] * df["sample_size"]

        logging.info("Filtering entries...")

        def apply_filter(df: pd.DataFrame, query_string: str) -> pd.DataFrame:
            included_indexes = df.query(query_string).index
            df["is_whitelisted_entry"] = df.index.isin(included_indexes)
            logging.debug(f"Whitelisting {len(included_indexes)}/{len(df)} rows where: {query_string}")
            return df

        df = apply_filter(df, f"allele_frequency > {HlaAlleleFreqNormaliser.NOT_PRESENT_FREQUENCY}")
        df = apply_filter(df, f"approx_observations > {HlaAlleleFreqNormaliser.MIN_ALLELE_OBSERVATIONS}")

        self._frequencies = df

    def normalise_frequencies(self) -> None:

        logging.info("Calculating normalised frequencies")

        groupings = self._frequencies.query("is_whitelisted_entry").groupby("allele")

        alleles = groupings.groups.keys()
        stats = pd.DataFrame(dict(
            locus=pd.Series(alleles).replace("\*.+$", "", regex=True).values,
            allele=alleles,
        ))

        stats["freq_mean_study_weighted"] = groupings.apply(
            lambda df: np.average(
                df['allele_frequency'],
                weights=df['study_weight'])
        ).values

        freq_mean_study_weighted_locus_sum = stats.groupby("locus")["freq_mean_study_weighted"].sum()
        stats["freq_mean_study_weighted_locus_sum"] = freq_mean_study_weighted_locus_sum[stats["locus"]].values

        stats["freq_normalised"] = stats["freq_mean_study_weighted"] / stats["freq_mean_study_weighted_locus_sum"]

        self._stats = stats

    @staticmethod
    def write_file(df: pd.DataFrame, output_path: str, *args, **kwargs) -> None:
        logging.info(f"Writing file: {output_path}")
        df.to_csv(output_path, index=False, *args, **kwargs)

    def write_debug_files(self) -> None:
        self.write_file(self._frequencies, f"{self.output_dir}/afnd.parsed.tsv.gz", sep="\t")
        self.write_file(self._stats, f"{self.output_dir}/afnd.stats.tsv.gz", sep="\t")

    def write_lilac_frequencies(self) -> None:

        logging.info("Writing debug files")

        lilac_frequencies = pd.DataFrame(dict(
            Allele = self._stats["allele"].values,
            CohortFrequency = self._stats["freq_normalised"].values,
        ))

        self.write_file(lilac_frequencies, f"{self.output_dir}/lilac_allele_frequencies.csv")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--input_file", help="Path to the TSV file containing the scraped AFND data")
    parser.add_argument("--output_dir", help="Output directory")
    parser.add_argument("--write_debug_files", help="Write tables containing intermediate calculations",
                        default=False, action=argparse.BooleanOptionalAction)

    args = parser.parse_args()

    input_file = args.input_file
    output_dir = args.output_dir
    write_debug_files = args.write_debug_files

    normaliser = HlaAlleleFreqNormaliser(afnd_path=input_file, output_dir=output_dir)

    normaliser.load_afnd()
    normaliser.normalise_frequencies()
    normaliser.write_lilac_frequencies()

    if write_debug_files:
        normaliser.write_debug_files()

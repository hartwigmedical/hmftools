"""
Adapted code from:
https://github.com/slowkow/allelefrequencies/blob/main/allelefrequencies.py
"""

from __future__ import annotations

import argparse
import logging

import os
import re
import gzip
import urllib
from enum import Enum

import pandas.errors
import requests
from bs4 import BeautifulSoup

import pandas as pd


logging.basicConfig(level=logging.DEBUG, format="%(asctime)s [%(levelname)-5s] %(message)s", datefmt="%H:%M:%S")
logging.addLevelName(logging.WARNING,  "WARN ")


class PopulationStandard(Enum):
    GOLD = "g"
    GOLD_AND_SILVER = "s"
    ALL = ""

    def get_file_string(self) -> str:
        if self == PopulationStandard.GOLD:
            return "gold_standard"

        if self == PopulationStandard.GOLD_AND_SILVER:
            return "gold_and_silver_standard"

        if self == PopulationStandard.ALL:
            return "any_standard"


class AlleleGroup(Enum):

    DIGITS_2 = 2
    DIGITS_4 = 4


class HlaLocusType(Enum):
    CLASSICAL = "Classical"
    NON_CLASSICAL = "Non-classical"

    def get_file_string(self) -> str:
        if self == HlaLocusType.CLASSICAL:
            return "classical"

        if self == HlaLocusType.NON_CLASSICAL:
            return "non_classical"


class DownloadableHlaLocus(Enum):

    A = "A"
    B = "B"
    C = "C"

    DPA1 = "DPA1"
    DPB1 = "DPB1"
    DQA1 = "DQA1"
    DQB1 = "DQB1"
    DRB1 = "DRB1"


def form_file_prefix(
        locus_or_string: DownloadableHlaLocus | str,
        hla_level: AlleleGroup,
        hla_locus_type: HlaLocusType,
        population_standard: PopulationStandard
) -> str:

    if isinstance(locus_or_string, DownloadableHlaLocus):
        locus_string = locus_or_string.value
    else:
        locus_string = locus_or_string

    prefix = f"{locus_string}.digits_{hla_level.value}.{hla_locus_type.get_file_string()}.{population_standard.get_file_string()}"
    return prefix


class HlaAlleleFreqDownloader:

    def __init__(
            self,
            output_dir: str,
            locus: DownloadableHlaLocus,
            hla_level: AlleleGroup = AlleleGroup.DIGITS_2,
            hla_locus_type: HlaLocusType = HlaLocusType.CLASSICAL,
            population_standard: PopulationStandard = PopulationStandard.GOLD
    ):
        self.locus = locus
        self.hla_level = hla_level
        self.hla_locus_type = hla_locus_type
        self.population_standard = population_standard

        if not os.path.exists(output_dir):
            logging.error(f"Invalid path: {output_dir}")
            raise Exception

        self.output_dir = output_dir
        self.cache_dir = f"{output_dir}/cache"

        self.main_url = self._form_main_url()

    def _form_main_url(self) -> str:
        base = 'http://www.allelefrequencies.net/hla6006a.asp'

        params = {
            'hla_locus': self.locus.value,
            'hla_level': self.hla_level.value,
            'hla_locus_type': self.hla_locus_type.value,
            'standard': self.population_standard.value
        }

        url = f"{base}?{urllib.parse.urlencode(params)}"

        return url

    @staticmethod
    def get_url(url: str, cache_dir: str) -> str:

        if not os.path.exists(cache_dir):
            logging.info(f"Making cache dir: {cache_dir}")
            os.makedirs(cache_dir, exist_ok=True)

        cache_basename = re.sub(r"[/\\?%*:|\"<>\x7F\x00-\x1F ]", "-", url)
        cache_file = f'{cache_dir}/{cache_basename}.html.gz'

        if not os.path.exists(cache_file):
            req = requests.get(url)
            text = req.text
            with gzip.open(cache_file, 'wt') as out:
                out.write(text)

        else:
            text = gzip.open(cache_file, 'rt').read()

        return text

    @staticmethod
    def get_all_pages_as_text(url: str, cache_dir: str) -> list[str]:

        logging.info(f'Scraping URL: {url}')

        text = HlaAlleleFreqDownloader.get_url(url, cache_dir)
        bs = BeautifulSoup(text, 'html.parser')

        # Number of pages of results
        n_pages = bs(text=re.compile(r' of \d'))
        if n_pages:
            n_pages = int(n_pages[0].strip()[3:].replace(',', ''))
        else:
            n_pages = 1

        if bs(text=re.compile(r'we did not find any results')):
            return [""]

        logging.info(f'Downloading content from {n_pages} pages')

        pages = []
        for i in range(1, n_pages + 1):
            text = HlaAlleleFreqDownloader.get_url(f'{url}&page={i}', cache_dir)
            pages.append(text)

        return pages

    @staticmethod
    def parse_one_page(text: str) -> pd.DataFrame:

        bs = BeautifulSoup(text, 'html.parser')

        table = bs.find('table', {'class': 'tblNormal'})
        trs = table.find_all('tr')

        rows = []
        for tr in trs[1:]:
            fields = [td.get_text(strip=True) for td in tr.find_all('td')]
            rows.append(fields)
        df = pd.DataFrame(rows)

        columns = dict(
            allele=1,
            population=3,
            perc_with_allele=4,
            allele_frequency=5,
            sample_size=7
        )

        df_subset = df[list(columns.values())]
        df_subset.columns = columns.keys()

        return df_subset

    def parse_all_pages(self, pages: list[str]) -> pd.DataFrame:

        logging.info("Combining pages to a single data frame")
        df_list = []

        for i in range(0, len(pages)):

            page_num = i+1
            logging.debug(f"Parsing page: {page_num} / {len(pages)}")

            text = pages[i]

            if len(text) == 0:
                logging.warning(f"Page {page_num} / {len(pages)} has no text to parse")
                continue

            df = HlaAlleleFreqDownloader.parse_one_page(text)
            df_list.append(df)

        if len(df_list) == 0:
            logging.warning(f"No content downloaded from URL: {self.main_url}")
            return pd.DataFrame()

        df_merged = pd.concat(df_list)
        return df_merged

    def download(self) -> pd.DataFrame:

        output_prefix = form_file_prefix(
            locus_or_string=self.locus,
            hla_level=self.hla_level,
            hla_locus_type=self.hla_locus_type,
            population_standard=self.population_standard
        )

        output_path = f"{self.output_dir}/{output_prefix}.tsv"

        if os.path.exists(output_path):
            logging.info(f"Loading cache file: {output_path}")

            try:
                return pd.read_csv(output_path, sep="\t")
            except pandas.errors.EmptyDataError:
                logging.warning(f"Skipping loading empty cache file: {output_path}")
                return pd.DataFrame()

        pages = self.get_all_pages_as_text(self.main_url, self.cache_dir)
        df = self.parse_all_pages(pages)

        if len(df) != 0:
            logging.info(f"Writing file: {output_path}")
        else:
            logging.info(f"Writing empty file: {output_path}")

        df.to_csv(output_path, index=False, sep="\t")

        return df


def download_all_loci(
        output_dir: str,
        loci: list[DownloadableHlaLocus],
        hla_level: AlleleGroup = AlleleGroup.DIGITS_2,
        hla_locus_type: HlaLocusType = HlaLocusType.CLASSICAL,
        population_standard: PopulationStandard = PopulationStandard.GOLD,
) -> pd.DataFrame:

    output_prefix = form_file_prefix(
        locus_or_string="afnd",
        hla_level=hla_level,
        hla_locus_type=hla_locus_type,
        population_standard=population_standard
    )

    output_path = f"{output_dir}/{output_prefix}.tsv.gz"

    df_list = []
    for locus in loci:
        logging.info(f"Downloading frequencies for: {locus.value}")

        output_subdir = f"{output_dir}/{locus.value}"
        os.makedirs(output_subdir, exist_ok=True)

        self = HlaAlleleFreqDownloader(
            output_dir = output_subdir,
            locus = locus,
            hla_level = hla_level,
            hla_locus_type=hla_locus_type,
            population_standard = population_standard
        )
        df = self.download()
        df_list.append(df)

        logging.info(" ")

    df = pd.concat(df_list)

    logging.info(f"Writing file: {output_path}")
    df.to_csv(output_path, sep="\t", index=False)

    return df


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--output_dir", help="Output directory")

    parser.add_argument(
        "--locus", default=None,
        help="Optional. Comma separated values of HLA loci to download (all if unspecified). "
             "Valid values: " + ",".join([enum.value for enum in DownloadableHlaLocus])
    )

    parser.add_argument(
        "--hla_level", default=AlleleGroup.DIGITS_2,
        help="Optional. HLA locus type. Valid values: " + ",".join([str(enum.value) for enum in AlleleGroup])
    )

    parser.add_argument(
        "--hla_locus_type", default=HlaLocusType.CLASSICAL,
        help="Optional. HLA locus type. Valid values: " + ",".join([enum.value for enum in HlaLocusType])
    )

    parser.add_argument(
        "--population_standard", default=PopulationStandard.GOLD,
        help="Optional. Population standard. Valid values: " + ",".join([enum.value for enum in PopulationStandard])
    )


    args = parser.parse_args()

    output_dir = args.output_dir
    hla_level = AlleleGroup(int(args.hla_level))
    hla_locus_type = HlaLocusType(args.hla_locus_type)
    population_standard = PopulationStandard(args.population_standard)

    if args.locus is not None:
        loci_strings = args.locus.split(",")
        loci = [DownloadableHlaLocus(string) for string in loci_strings]
    else:
        loci = DownloadableHlaLocus

    download_all_loci(
        output_dir=output_dir,
        loci=loci,
        hla_level=hla_level,
        hla_locus_type=hla_locus_type,
        population_standard=population_standard,
    )


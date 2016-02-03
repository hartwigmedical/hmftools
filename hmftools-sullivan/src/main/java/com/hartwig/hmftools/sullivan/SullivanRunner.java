package com.hartwig.hmftools.sullivan;

public class SullivanRunner {

    public static void main(String[] args) {
        String originalFastqPath =
                "/Users/kduyvesteyn/hmf/data/100k_reads_hiseq/Sample_CONTROLP25_H7YRLADXX/CONTROLP25_H7YRLADXX_ATCACG_L001_R1_001.subset_100k.fastq";
        String recreatedFastqPath =
                "/Users/kduyvesteyn/hmf/out/fastq/CONTROLP25_H7YRLADXX_ATCACG_L001_1.fastq";

        SullivanAlgo.runSullivanAlgo(originalFastqPath, recreatedFastqPath);
    }
}

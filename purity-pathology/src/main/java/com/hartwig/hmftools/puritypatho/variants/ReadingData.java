package com.hartwig.hmftools.puritypatho.variants;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ReadingData {
    private static final Logger LOGGER = LogManager.getLogger(ReadingData.class);

    private static final String EXTENSION = ".amber.baf";
    private static final String AMBER_DIR = "amber";
    private static final String MAP_BAF_FILE = "/data/common/dbs/CytoScanHD/CytoScanHD_hg19_SNPs_sorted.bed";

    public static void readingFiles(@NotNull String runsFolderPath, @NotNull String tumorSample, @NotNull String countSet) throws IOException {
        LOGGER.info("Creating output file");
        String fileName = VariantDetection.generateOutputFile();

        LOGGER.info("Reading CytoScan file: " + MAP_BAF_FILE);
        final List<String> readingCytoScanFile = CytoScanFile.read(MAP_BAF_FILE);
        Multimap<String, String> multimapCyto = VariantDetection.extractCytoData(readingCytoScanFile);

        LOGGER.info("Reading Amber file: " + runsFolderPath + File.separator + AMBER_DIR + File.separator + tumorSample + EXTENSION);
        final String amberFile = runsFolderPath + File.separator + AMBER_DIR + File.separator + tumorSample + EXTENSION;
        final Multimap<String, AmberBAF> amberBAFs = AmberBAFFile.read(amberFile);
        List<AmberBAF> sortedBafs = Lists.newArrayList(amberBAFs.values());
        Collections.sort(sortedBafs);

        LOGGER.info("Total BAF points of file from sample " + tumorSample + ": " + sortedBafs.size());
        final List<String> finalPurityData = AmberBAFFile.readingPurityData(sortedBafs);
        VariantDetection.extractAmberData(finalPurityData, multimapCyto, fileName, countSet);
    }
}
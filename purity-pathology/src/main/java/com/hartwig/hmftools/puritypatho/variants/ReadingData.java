package com.hartwig.hmftools.puritypatho.variants;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.jetbrains.annotations.NotNull;

public class ReadingData {
    private static final String EXTENSION = ".amber.baf";
    private static final String AMBER_DIR = "amber";
    private static final String MAP_BAF_FILE = "/data/common/dbs/CytoScanHD/CytoScanHD_hg19_SNPs_sorted.bed";
    private static final Logger LOGGER = LogManager.getLogger(ReadingData.class);

    public static void readingSetAmber(@NotNull String runsFolderPath, @NotNull String tumorSample) throws IOException{
        LOGGER.info("reading Amber file: " + runsFolderPath + File.separator + AMBER_DIR + File.separator + tumorSample + EXTENSION);
        final String Amberfile = runsFolderPath + File.separator + AMBER_DIR + File.separator + tumorSample + EXTENSION;
        final Multimap<String, AmberBAF> AmberfileReading = AmberBAFFile.read(Amberfile);
        VariantDetection.checkVariants(AmberfileReading);
    }

    public static void readingCyto () throws IOException {
        LOGGER.info("reading Cytoscan file" + MAP_BAF_FILE);
        final List CytoFileReading = Files.readAllLines(new File(MAP_BAF_FILE).toPath());
        VariantDetection.runCyto(CytoFileReading);
    }
}

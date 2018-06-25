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

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public class ReadingData {
    public static final String DELIMITER = "\t";
    public static final String EXTENSION = ".amber.baf";
    public static final String AMBER_DIR = "amber";
    public static final String MAP_BAF_FILE = "/data/common/dbs/CytoScanHD/CytoScanHD_hg19_SNPs_sorted.bed";
    public static final Logger LOGGER = LogManager.getLogger(ReadingData.class);

    public static void readingFiles(@NotNull String runsFolderPath, @NotNull String tumorSample) throws IOException {
        final String Amberfile = runsFolderPath + File.separator + AMBER_DIR + File.separator + tumorSample + EXTENSION;

        final Multimap<String, AmberBAF> AmberBAFFiles = AmberBAFFile.read(Amberfile);
        List<AmberBAF> sortedBafs = Lists.newArrayList(AmberBAFFiles.values());
        Collections.sort(sortedBafs);

        final List<String> finalPurityData = AmberBAFFile.readingPurityData(sortedBafs);
        LOGGER.info("reading Amber file: " + runsFolderPath + File.separator + AMBER_DIR + File.separator + tumorSample + EXTENSION);

        for (String lineAmber : finalPurityData) {
            String[] partsAmber = lineAmber.split(DELIMITER);
            String chromosomesAmber = partsAmber[0];
            String positionsAmber = partsAmber[1];
            LOGGER.info(lineAmber);
        }

        LOGGER.info("MAP_BAF_FILE: " + MAP_BAF_FILE);
        final List<String> readingCytoScanFile = CytoScanFile.read(MAP_BAF_FILE);
        for (String lineCyto : readingCytoScanFile) {
            String[] partsCyto = lineCyto.split(DELIMITER);
            String chromosomesCyto = partsCyto[0];
            String positionsStartCyto = partsCyto[1];
            String positionsEndCyto = partsCyto[2];
            LOGGER.info(lineCyto);
        }
    }
}
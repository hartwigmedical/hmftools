package com.hartwig.hmftools.puritypatho.variants;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.ListMultimap;
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
    private static final String EXTENSION = ".amber.baf";
    private static final String AMBER_DIR = "amber";
    private static final String MAP_BAF_FILE = "/data/common/dbs/CytoScanHD/CytoScanHD_hg19_SNPs_sorted.bed";
    private static final Logger LOGGER = LogManager.getLogger(ReadingData.class);

    public static void readingFiles(@NotNull String runsFolderPath, @NotNull String tumorSample) throws IOException {
        LOGGER.info("Creating output file");
        String fileName = VariantDetection.GenerateOutpurFile();

        LOGGER.info("Reading CytoScan file");
        final List<String> readingCytoScanFile = CytoScanFile.read(MAP_BAF_FILE);
        ListMultimap<String, String> multimapCyto = VariantDetection.ExtractCytoData(readingCytoScanFile);

        LOGGER.info("Reading Amber file: " + runsFolderPath + File.separator + AMBER_DIR + File.separator + tumorSample + EXTENSION);
        final String Amberfile = runsFolderPath + File.separator + AMBER_DIR + File.separator + tumorSample + EXTENSION;
        final Multimap<String, AmberBAF> AmberBAFFiles = AmberBAFFile.read(Amberfile);
        List<AmberBAF> sortedBafs = Lists.newArrayList(AmberBAFFiles.values());
        Collections.sort(sortedBafs);

        LOGGER.info("Total BAF points of file from sample " + tumorSample + ": " + sortedBafs.size());
        final List<String> finalPurityData = AmberBAFFile.readingPurityData(sortedBafs);
        VariantDetection.ExtractAmberData(finalPurityData, multimapCyto, fileName);
    }
}
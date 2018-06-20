package com.hartwig.hmftools.puritypatho.variants;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.jetbrains.annotations.NotNull;

public class ReadingData {
    private static final String EXTENSION = ".amber.baf";
    private static final String AMBER_DIR = "amber";
    private static final String MAP_BAF_FILE = "/data/common/dbs/CytoScanHD/CytoScanHD_hg19_SNPs_sorted.bed";
    private static final Logger LOGGER = LogManager.getLogger(ReadingData.class);

    public static void readingFiles(@NotNull String runsFolderPath, @NotNull String tumorSample) throws IOException{
        LOGGER.info("reading Amber file: " + runsFolderPath + File.separator + AMBER_DIR + File.separator + tumorSample + EXTENSION);
        final String Amberfile = runsFolderPath + File.separator + AMBER_DIR + File.separator + tumorSample + EXTENSION;
        int countAmber = readingLines(Amberfile).split("\n").length;

        LOGGER.info("reading Cytoscan file" + MAP_BAF_FILE);
        int countCyto = readingLines(MAP_BAF_FILE).split("\n").length;

        VariantDetection.checkVariants(readingLines(Amberfile), countAmber,
                readingLines(MAP_BAF_FILE), countCyto);
    }

    public static String readingLines (@NotNull String file) throws IOException {
        FileReader fileReader = new FileReader(file);
        BufferedReader bufferedReader = new BufferedReader(fileReader);
        StringBuffer stringBuffer = new StringBuffer();
        String line;
        while ((line = bufferedReader.readLine()) != null) {
            stringBuffer.append(line);
            stringBuffer.append("\n");
        }
        fileReader.close();
        return stringBuffer.toString();
    }
}

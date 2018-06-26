package com.hartwig.hmftools.puritypatho.variants;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
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
    public static final String DELIMITER = "\t";
    public static final String EXTENSION = ".amber.baf";
    public static final String AMBER_DIR = "amber";
    public static final String MAP_BAF_FILE = "/data/common/dbs/CytoScanHD/CytoScanHD_hg19_SNPs_sorted.bed";
    public static final Logger LOGGER = LogManager.getLogger(ReadingData.class);

    public static void filterVariant(@NotNull String chromosome, @NotNull String chromosomesAmber, @NotNull String positionsAmber,
            @NotNull ListMultimap<String, String> multimapCyto) {
        LOGGER.info(chromosome + ": "  + chromosomesAmber + "\t" + positionsAmber);
        LOGGER.info(multimapCyto.get(chromosomesAmber).contains(positionsAmber));
    }

    public static void readingFiles(@NotNull String runsFolderPath, @NotNull String tumorSample) throws IOException {

        LOGGER.info("MAP_BAF_FILE: " + MAP_BAF_FILE);
        ListMultimap<String, String> multimapCyto = ArrayListMultimap.create();
        final List<String> readingCytoScanFile = CytoScanFile.read(MAP_BAF_FILE);
        for (String lineCyto : readingCytoScanFile) {
            String[] partsCyto = lineCyto.split(DELIMITER);
            String chromosomesCyto = partsCyto[0];
            String positionsCyto = partsCyto[1];
            multimapCyto.put(chromosomesCyto, positionsCyto);
        }

        LOGGER.info("reading Amber file: " + runsFolderPath + File.separator + AMBER_DIR + File.separator + tumorSample + EXTENSION);
        final String Amberfile = runsFolderPath + File.separator + AMBER_DIR + File.separator + tumorSample + EXTENSION;
        final Multimap<String, AmberBAF> AmberBAFFiles = AmberBAFFile.read(Amberfile);
        List<AmberBAF> sortedBafs = Lists.newArrayList(AmberBAFFiles.values());
        Collections.sort(sortedBafs);

        final List<String> finalPurityData = AmberBAFFile.readingPurityData(sortedBafs);

        for (String lineAmber : finalPurityData) {
            LOGGER.info("lineAmber: " + lineAmber);
            String[] partsAmber = lineAmber.split(DELIMITER);
            String chromosomesAmber = partsAmber[0];
            String positionsAmber = partsAmber[1];

           if (chromosomesAmber.equals("1")) {
               filterVariant("1", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("2")) {
               filterVariant("2", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("3")) {
               filterVariant("3", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("4")) {
               filterVariant("4", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("5")) {
               filterVariant("5", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("6")) {
               filterVariant("6", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("7")) {
               filterVariant("7", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("8")) {
               filterVariant("8", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("9")) {
               filterVariant("9", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("10")) {
               filterVariant("10", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("11")) {
                filterVariant("11", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("12")) {
                filterVariant("12", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("13")) {
                filterVariant("13", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("14")) {
                filterVariant("14", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("15")) {
                filterVariant("15", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("16")) {
                filterVariant("16", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("17")) {
                filterVariant("17", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("18")) {
                filterVariant("18", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("19")) {
                filterVariant("19", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("20")) {
                filterVariant("20", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("21")) {
                filterVariant("21", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("22")) {
                filterVariant("22", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("X")) {
                filterVariant("X", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("Y")) {
                filterVariant("Y", chromosomesAmber, positionsAmber, multimapCyto);
           } else if (chromosomesAmber.equals("MT")) {
                filterVariant("MT", chromosomesAmber, positionsAmber, multimapCyto);
           } else {
                LOGGER.info("No known chromosome value!");
           }
        }
    }
}
package com.hartwig.hmftools.puritypatho.variants;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public class VariantDetection {
    private static final Logger LOGGER = LogManager.getLogger(VariantDetection.class);
    private static final String DELIMITER = "\t";

    public static String GenerateOutpurFile() throws IOException{
        final String filename = WritingData.generateOutputFileName();
        WritingData.writeToFileHeader(filename);
        return filename;
    }

    public static ListMultimap<String, String> ExtractCytoData(@NotNull List<String> readingCytoScanFile) {
        ListMultimap<String, String> multimapCyto = ArrayListMultimap.create();
        for (String lineCyto : readingCytoScanFile) {
            String[] partsCyto = lineCyto.split(DELIMITER);
            String chromosomesCyto = partsCyto[0];
            String positionsCyto = partsCyto[1];
            multimapCyto.put(chromosomesCyto, positionsCyto);
        }
        return multimapCyto;
    }

    public static void ExtractAmberData(@NotNull List<String> finalPurityData, @NotNull ListMultimap<String, String> multimapCyto,
            @NotNull String fileName) throws IOException {
        for (String lineAmber : finalPurityData) {
            String[] partsAmber = lineAmber.split(DELIMITER);
            String chromosomesAmber = partsAmber[0];
            String positionsAmber = partsAmber[1];
            int countAmber = 0;

            if (chromosomesAmber.equals("1")) {
                filterVariant("1", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("2")) {
                filterVariant("2", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("3")) {
                filterVariant("3", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("4")) {
                filterVariant("4", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("5")) {
                filterVariant("5", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("6")) {
                filterVariant("6", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("7")) {
                filterVariant("7", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("8")) {
                filterVariant("8", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("9")) {
                filterVariant("9", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("10")) {
                filterVariant("10", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("11")) {
                filterVariant("11", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("12")) {
                filterVariant("12", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("13")) {
                filterVariant("13", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("14")) {
                filterVariant("14", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("15")) {
                filterVariant("15", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("16")) {
                filterVariant("16", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("17")) {
                filterVariant("17", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("18")) {
                filterVariant("18", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("19")) {
                filterVariant("19", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("20")) {
                filterVariant("20", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("21")) {
                filterVariant("21", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("22")) {
                filterVariant("22", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("X")) {
                filterVariant("X", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("Y")) {
                filterVariant("Y", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else if (chromosomesAmber.equals("MT")) {
                filterVariant("MT", chromosomesAmber, positionsAmber, multimapCyto, lineAmber, countAmber, fileName);
            } else {
                LOGGER.info("No known chromosome value!");
            }
        }
    }

    public static void filterVariant(@NotNull String chromosome, @NotNull String chromosomesAmber, @NotNull String positionsAmber,
            @NotNull ListMultimap<String, String> multimapCyto, @NotNull String lineAmber, int countAmber, String fileName) throws
            IOException {
        if (multimapCyto.get(chromosomesAmber).contains(positionsAmber)) {
            countAmber ++;
            LOGGER.info(chromosomesAmber + DELIMITER + positionsAmber + DELIMITER + countAmber);
            WritingData.writeToFile(fileName, chromosomesAmber , positionsAmber, countAmber);
        } else {
            WritingData.writeToFile(fileName, chromosomesAmber , positionsAmber, countAmber);
        }
    }

}

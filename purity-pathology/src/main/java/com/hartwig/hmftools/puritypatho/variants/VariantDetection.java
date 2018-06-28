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
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("2")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("3")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("4")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("5")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("6")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("7")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("8")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("9")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("10")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("11")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("12")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("13")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("14")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("15")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("16")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("17")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("18")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("19")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("20")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("21")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("22")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("X")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("Y")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else if (chromosomesAmber.equals("MT")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName);
            } else {
                LOGGER.info("No known chromosome value!");
            }
        }
    }

    public static void filterVariant(@NotNull String chromosomesAmber, @NotNull String positionsAmber,
            @NotNull ListMultimap<String, String> multimapCyto,int countAmber, String fileName) throws
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

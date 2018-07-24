package com.hartwig.hmftools.puritypatho.variants;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class VariantDetection {
    private static final Logger LOGGER = LogManager.getLogger(VariantDetection.class);
    private static final String DELIMITER = "\t";

    @NotNull
    public static String generateOutputFile() {
        return WritingData.generateOutputFileName();
    }

    @NotNull
    public static Multimap<String, String> extractCytoData(@NotNull List<String> readingCytoScanFile) {
        Multimap<String, String> multimapCyto = ArrayListMultimap.create();
        for (String lineCyto : readingCytoScanFile) {
            String[] partsCyto = lineCyto.split(DELIMITER);
            String chromosomesCyto = partsCyto[0];
            String positionsCyto = partsCyto[2];
            multimapCyto.put(chromosomesCyto, positionsCyto);
        }
        return multimapCyto;
    }

    @NotNull
    private static Map<String, String> readingOutput(@NotNull String fileName) throws IOException {
        final List<String> output = ReadingFileVariantDetection.read(fileName);
        WritingData.writeToFileHeader(fileName);
        final Map<String, String> resultOutput = Maps.newHashMap();

        for (String lineOutput : output) {
            String[] partsOutput = lineOutput.split(DELIMITER);
            String outputGenomic = partsOutput[0] + "," + partsOutput[1];
            String outputCount = partsOutput[2];
            resultOutput.put(outputGenomic, outputCount);
        }
        return resultOutput;
    }

    public static void extractAmberData(@NotNull Multimap<String, String> multimapCyto,
            @NotNull String fileName, @NotNull String countSet, List<AmberBAF> amberBAFS) throws IOException {
        if (countSet.equals("1")) {
            WritingData.writeToFileHeader(fileName);
        }
        Map<String, String> resultOutput = readingOutput(fileName);
        Set<String> genomicPosition = resultOutput.keySet();

        for (AmberBAF baf : amberBAFS) {
            String chromosomesAmber = baf.chromosome();
            String positionsAmber = Long.toString(baf.position());
            int countAmber = 0;

            if (validChromosome(chromosomesAmber)) {
                filterVariant(chromosomesAmber,
                        positionsAmber,
                        multimapCyto,
                        countAmber,
                        fileName,
                        resultOutput,
                        genomicPosition,
                        countSet);
            } else {
                LOGGER.info("No known chromosome value!");
            }
        }
        uniqueValuesOfPreviousFile(genomicPosition, countSet, resultOutput, fileName);
    }

    private static Boolean validChromosome (@NotNull String chromosomesAmber) {
        List<String> chromosomesValids = Arrays.asList("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
                "18", "19", "20", "21", "22", "X", "Y", "MT");
       return chromosomesValids.contains(chromosomesAmber);
    }

    private static void uniqueValuesOfPreviousFile(@NotNull Set<String> genomicPosition, @NotNull String countSet,
            @NotNull Map<String, String> resultOutput, @NotNull String fileName) throws IOException {
        genomicPosition.remove("chromosome" + "," + "position");
        List<String> sortedGenomicPosition = Lists.newArrayList(genomicPosition);
        Collections.sort(sortedGenomicPosition);
        if (!countSet.equals("1")) {
            for (String position : sortedGenomicPosition) {
                String[] outputGenomic = position.split(",");
                Integer countValue = Integer.valueOf(resultOutput.get(position));
                WritingData.writeToFile(fileName, outputGenomic[0], outputGenomic[1], countValue);
            }
        }
    }

    private static void filterVariant(@NotNull String chromosomesAmber, @NotNull String positionsAmber,
            @NotNull Multimap<String, String> multimapCyto, int countAmber, @NotNull String fileName,
            @NotNull Map<String, String> resultOutput, @NotNull Set genomicPosition, @NotNull String countSet) throws IOException {
        if (multimapCyto.get(chromosomesAmber).contains(positionsAmber)) {
            if (countSet.equals("1")) {
                countAmber++;
                WritingData.writeToFile(fileName, chromosomesAmber, positionsAmber, countAmber);
            } else {
                final String position = chromosomesAmber + "," + positionsAmber;
                final boolean foundGenomicPositionInFile = genomicPosition.contains(position);
                if (foundGenomicPositionInFile) {
                    final Integer valueCount = Integer.valueOf(resultOutput.get(position));
                    int countCombined = valueCount + 1;
                    WritingData.writeToFile(fileName, chromosomesAmber, positionsAmber, countCombined);
                    genomicPosition.remove(chromosomesAmber + "," + positionsAmber);
                } else {
                    countAmber++;
                    WritingData.writeToFile(fileName, chromosomesAmber, positionsAmber, countAmber);
                }
            }
        }
    }
}

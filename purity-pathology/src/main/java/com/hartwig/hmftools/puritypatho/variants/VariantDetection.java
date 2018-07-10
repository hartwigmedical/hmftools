package com.hartwig.hmftools.puritypatho.variants;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;

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
    private static Multimap<String, String> readingOutput(@NotNull String fileName) throws IOException {
        final List<String> output = ReadingFileVariantDetection.read(fileName);
        WritingData.writeToFileHeader(fileName);
        final Multimap<String, String> resultOutput = ArrayListMultimap.create();

        for (String lineOutput : output) {
            String[] partsOutput = lineOutput.split(DELIMITER);
            String outputGenomic = partsOutput[0] + "," + partsOutput[1];
            String outputCount = partsOutput[2];
            resultOutput.put(outputGenomic, outputCount);
        }
        return resultOutput;
    }

    public static void extractAmberData(@NotNull List<String> finalPurityData, @NotNull Multimap<String, String> multimapCyto,
            @NotNull String fileName, @NotNull String countSet) throws IOException {
        if (countSet.equals("1")) {
            WritingData.writeToFileHeader(fileName);
        }
        Multimap<String, String> resultOutput = readingOutput(fileName);
        Set<String> genomicPosition = resultOutput.keySet();

        for (String lineAmber : finalPurityData) {
            String[] partsAmber = lineAmber.split(DELIMITER);
            String chromosomesAmber = partsAmber[0];
            String positionsAmber = partsAmber[1];
            int countAmber = 0;

            switch (chromosomesAmber) {
                case "1":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "2":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "3":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "4":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "5":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "6":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "7":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "8":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "9":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "10":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "11":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "12":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "13":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "14":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "15":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "16":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "17":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "18":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "19":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "20":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "21":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "22":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "X":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "Y":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                case "MT":
                    filterVariant(chromosomesAmber,
                            positionsAmber,
                            multimapCyto,
                            countAmber,
                            fileName,
                            resultOutput,
                            genomicPosition,
                            countSet);
                    break;
                default:
                    LOGGER.info("No known chromosome value!");
                    break;
            }
        }
        uniqueValuesOfPreviousFile(genomicPosition, countSet, resultOutput, fileName);
    }

    private static void uniqueValuesOfPreviousFile(@NotNull Set<String> genomicPosition, @NotNull String countSet,
            @NotNull Multimap<String, String> resultOutput, @NotNull String fileName) throws IOException {
        genomicPosition.remove("chromosome" + "," + "position");
        List<String> sortedGenomicPosition = Lists.newArrayList(genomicPosition);
        Collections.sort(sortedGenomicPosition);
        if (!countSet.equals("1")) {
            for (String position : sortedGenomicPosition) {
                String[] outputGenomic = position.split(",");
                String countValue = resultOutput.get(position).toString().replace("[", "");
                Integer countValueDef = Integer.valueOf(countValue.replace("]", ""));
                WritingData.writeToFile(fileName, outputGenomic[0], outputGenomic[1], countValueDef);
            }
        }
    }

    private static void filterVariant(@NotNull String chromosomesAmber, @NotNull String positionsAmber,
            @NotNull Multimap<String, String> multimapCyto, int countAmber, @NotNull String fileName,
            @NotNull Multimap<String, String> resultOutput, @NotNull Set genomicPosition, @NotNull String countSet) throws IOException {
        if (multimapCyto.get(chromosomesAmber).contains(positionsAmber)) {
            if (countSet.equals("1")) {
                countAmber++;
                WritingData.writeToFile(fileName, chromosomesAmber, positionsAmber, countAmber);
            } else {
                final String position = chromosomesAmber + "," + positionsAmber;
                final boolean foundGenomicPositionInFile = genomicPosition.contains(position);
                if (foundGenomicPositionInFile) {
                    final String valueCount = resultOutput.get(position).toString();
                    String valueCountNew = valueCount.replace("[", "");
                    int valueCountNewDef = Integer.valueOf(valueCountNew.replace("]", ""));
                    int countCombined = valueCountNewDef + 1;
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

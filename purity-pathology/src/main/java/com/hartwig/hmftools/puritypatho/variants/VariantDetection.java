package com.hartwig.hmftools.puritypatho.variants;

import java.io.IOException;

import java.util.List;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jfree.util.StringUtils;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public class VariantDetection {
    private static final Logger LOGGER = LogManager.getLogger(VariantDetection.class);
    private static final String DELIMITER = "\t";

    public static String GenerateOutputFile(@NotNull String countSet) throws IOException{
        final String filename = WritingData.generateOutputFileName();
        return filename;
    }

    public static ListMultimap<String, String> ExtractCytoData(@NotNull List<String> readingCytoScanFile) {
        ListMultimap<String, String> multimapCyto = ArrayListMultimap.create();
        for (String lineCyto : readingCytoScanFile) {
            String[] partsCyto = lineCyto.split(DELIMITER);
            String chromosomesCyto = partsCyto[0];
            String positionsCyto = partsCyto[2];
            multimapCyto.put(chromosomesCyto, positionsCyto);
        }
        return multimapCyto;
    }

    public static void ExtractAmberData(@NotNull List<String> finalPurityData, @NotNull ListMultimap<String, String> multimapCyto,
            @NotNull String fileName, @NotNull String countSet) throws IOException {
        if (countSet.equals("1")){
            WritingData.writeToFileHeader(fileName);
        }
        final List<String> output = ReadingFileVariantDetection.read(fileName);
        WritingData.writeToFileHeader(fileName);
        final ListMultimap<String, String> resultOutput = ArrayListMultimap.create();

        for (String lineOutput : output) {
            String[] partsOutput = lineOutput.split(DELIMITER);
            String outputGenomic = partsOutput[0] + "," + partsOutput[1];
            String outputCount = partsOutput[2];
            resultOutput.put(outputGenomic, outputCount);
        }
        Set genomicPosition = resultOutput.keySet();
        for (String lineAmber : finalPurityData) {
            String[] partsAmber = lineAmber.split(DELIMITER);
            String chromosomesAmber = partsAmber[0];
            String positionsAmber = partsAmber[1];
            int countAmber = 0;

            if (chromosomesAmber.equals("1")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("2")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("3")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("4")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("5")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("6")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("7")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("8")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("9")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("10")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("11")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("12")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("13")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("14")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("15")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("16")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("17")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("18")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("19")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("20")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("21")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("22")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("X")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("Y")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else if (chromosomesAmber.equals("MT")) {
                filterVariant(chromosomesAmber, positionsAmber, multimapCyto, countAmber, fileName, resultOutput, genomicPosition);
            } else {
                LOGGER.info("No known chromosome value!");
            }
        }
    }

    private static void filterVariant(@NotNull String chromosomesAmber, @NotNull String positionsAmber,
            @NotNull ListMultimap<String, String> multimapCyto,int countAmber, String fileName,
            @NotNull Multimap<String, String> resultOutput, @NotNull Set genomicPosition) throws
            IOException {
        if (multimapCyto.get(chromosomesAmber).contains(positionsAmber)) {
            final String position = chromosomesAmber + "," + positionsAmber;
            final Boolean foundGenomicPostionInFile = genomicPosition.contains(position);
            if (foundGenomicPostionInFile){
                final String value = resultOutput.get(position).toString();
                String valueNew = value.replace("[", "");
                int valueDef = Integer.valueOf(valueNew.replace("]", ""));
                int countCombined = valueDef + 1;
            }
            countAmber ++;
            WritingData.writeToFile(fileName, chromosomesAmber , positionsAmber, countAmber);
        }
    }
}

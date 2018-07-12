package com.hartwig.hmftools.puritypatho.variants;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class WritingData {
    private static final String pathName = "/home/lieke/analysis/";
    private static final String outputName = "variants.txt";
    private static final String DELIMITER = "\t";

    @NotNull
    public static String generateOutputFileName(){
        return pathName + outputName;
    }

    public static void writeToFile(@NotNull String fileName, @NotNull String chromosomesAmber, @NotNull String positionsAmber, int countAmber) throws
            IOException {
        Files.write(new File(fileName).toPath(), toLines(chromosomesAmber, positionsAmber, countAmber), StandardOpenOption.APPEND);
    }

    @NotNull
    private static List toLines(@NotNull String chromosomesAmber, @NotNull String positionsAmber, int countAmber) {
        final List<String> result = Lists.newArrayList();
        result.add(chromosomesAmber + DELIMITER + positionsAmber + DELIMITER + countAmber);
        return result;
    }

    public static void writeToFileHeader(@NotNull String fileName) throws IOException{
        Files.write(new File(fileName).toPath(), toLinesHeader());
    }

    @NotNull
    private static List toLinesHeader() {
        final List<String> result = Lists.newArrayList();
        result.add("chromosome" + DELIMITER + "position" + DELIMITER + "count");
        return result;
    }
}

package com.hartwig.hmftools.puritypatho.variants;

import com.google.common.collect.Lists;

import java.io.File;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

public class VariantDetection {
    private static final Logger LOGGER = LogManager.getLogger(VariantDetection.class);
    private static final String pathName = "/home/lieke/analysis/";
    private static final String outputName = "variants.txt";

    public static void checkVariants(@NotNull String AmberfileReading, @NotNull String CytoFileReading){
    }

    public static String generateOutputFileName(){
        return pathName + outputName;
    }

    public static void write (@NotNull final String fileName, @NotNull final String results) throws IOException {
        Files.write(new File(fileName).toPath(), toLines(results));
    }

    @NotNull
    private static List<String> toLines(@NotNull final String results) {
        final List<String> result = Lists.newArrayList();
        result.add("chromosome" + "position1" + "position2" + "counts");
        result.add(results);
        return result;
    }
}

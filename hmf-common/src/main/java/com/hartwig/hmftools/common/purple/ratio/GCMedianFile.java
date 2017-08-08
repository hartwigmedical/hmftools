package com.hartwig.hmftools.common.purple.ratio;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public enum GCMedianFile {
    ;

    private static final String DELIMITER = "\t";
    private static final String EXTENSION = ".purple.gc.median";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(@NotNull final String fileName, @NotNull final List<GCMedian> medians) throws IOException {
        Files.write(new File(fileName).toPath(), toLines(medians));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<GCMedian> medians) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        medians.stream().map(GCMedianFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("GCContent").add("MedianReadCount").toString();
    }

    @NotNull
    private static String toString(@NotNull final GCMedian median) {
        return new StringJoiner(DELIMITER).add(String.valueOf(median.gcContent())).add(String.valueOf(median.medianCount())).toString();
    }

}

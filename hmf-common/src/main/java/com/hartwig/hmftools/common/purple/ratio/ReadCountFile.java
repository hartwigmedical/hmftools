package com.hartwig.hmftools.common.purple.ratio;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class ReadCountFile {

    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "Chr";
    private static final String EXTENSION = ".cobalt.ratio";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void createFile(@NotNull final String filename) throws IOException {
        Files.write(new File(filename).toPath(), (header() + "\n"). getBytes());
    }

    public static void append(@NotNull final String filename, @NotNull final List<ReadCount> readCounts) throws IOException {
        Files.write(new File(filename).toPath(), toLines(readCounts), StandardOpenOption.APPEND);
    }

    @NotNull
    static List<String> toLines(@NotNull final List<ReadCount> readCounts) {
        final List<String> lines = Lists.newArrayList();
        readCounts.stream().map(ReadCountFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("Chromosome").add("Position").add("ReadCount").toString();
    }

    @NotNull
    private static String toString(@NotNull final ReadCount ratio) {
        return new StringJoiner(DELIMITER).add(String.valueOf(ratio.chromosome()))
                .add(String.valueOf(ratio.position()))
                .add(String.valueOf(ratio.readCount()))
                .toString();
    }
}

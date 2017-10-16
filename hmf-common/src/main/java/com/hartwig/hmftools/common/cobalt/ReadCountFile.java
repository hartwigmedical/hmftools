package com.hartwig.hmftools.common.cobalt;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;

import org.jetbrains.annotations.NotNull;

public class ReadCountFile {

    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "#";
    private static final String EXTENSION = ".cobalt.count";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @Deprecated
    public static void createFile(int windowSize, @NotNull final String filename) throws IOException {
        Files.write(new File(filename).toPath(), (metaData(windowSize) + "\n" + header() + "\n").getBytes());
    }

    @Deprecated
    public static void append(@NotNull final String filename, @NotNull final List<ReadCount> readCounts) throws IOException {
        Files.write(new File(filename).toPath(), toLines(readCounts), StandardOpenOption.APPEND);
    }

    @NotNull
    public static Multimap<String, ReadCount> readFile(@NotNull final String filename) throws IOException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<ReadCount> readCounts) {
        final List<String> lines = Lists.newArrayList();
        readCounts.stream().map(ReadCountFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String metaData(int windowSize) {
        return HEADER_PREFIX + "#windowSize=" + windowSize;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, HEADER_PREFIX, "").add("Chromosome").add("Position").add("ReadCount").toString();
    }

    @NotNull
    private static String toString(@NotNull final ReadCount ratio) {
        return new StringJoiner(DELIMITER).add(String.valueOf(ratio.chromosome()))
                .add(String.valueOf(ratio.position()))
                .add(String.valueOf(ratio.readCount()))
                .toString();
    }

    @NotNull
    private static Multimap<String, ReadCount> fromLines(@NotNull final List<String> lines) {

        final Multimap<String, ReadCount> result = ArrayListMultimap.create();
        for (String line : lines) {
            if (!line.startsWith(HEADER_PREFIX)) {
                final ReadCount readCount = fromLine(line);
                result.put(readCount.chromosome(), readCount);
            }
        }

        return result;
    }

    @NotNull
    private static ReadCount fromLine(@NotNull final String ratioLine) {
        final String[] values = ratioLine.split(DELIMITER);

        final String chromosome = values[0].trim();
        final long position = Long.valueOf(values[1].trim());
        final int readCount = Integer.valueOf(values[2].trim());

        return ImmutableReadCount.builder().chromosome(chromosome).position(position).readCount(readCount).build();
    }

}

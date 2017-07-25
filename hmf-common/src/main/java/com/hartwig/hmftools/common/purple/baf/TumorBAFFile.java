package com.hartwig.hmftools.common.purple.baf;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;

import org.jetbrains.annotations.NotNull;

public enum TumorBAFFile {
    ;
    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "#";
    private static final String EXTENSION = ".purple.baf";

    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static Multimap<String, TumorBAF> read(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull Multimap<String, TumorBAF> bafs) throws IOException {
        List<TumorBAF> sortedBafs = Lists.newArrayList(bafs.values());
        Collections.sort(sortedBafs);
        write(filename, sortedBafs);
    }

    public static void write(@NotNull final String filename, @NotNull List<TumorBAF> bafs) throws IOException {
        Files.write(new File(filename).toPath(), toLines(bafs));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<TumorBAF> purity) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        purity.stream().map(TumorBAFFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, HEADER_PREFIX, "").add("chromosome").add("position").add("baf").toString();
    }

    @NotNull
    private static String toString(@NotNull final TumorBAF ratio) {
        return new StringJoiner(DELIMITER).add(String.valueOf(ratio.chromosome()))
                .add(String.valueOf(ratio.position()))
                .add(String.valueOf(ratio.baf()))
                .toString();
    }

    @NotNull
    private static Multimap<String, TumorBAF> fromLines(@NotNull List<String> lines) {
        Multimap<String, TumorBAF> result = ArrayListMultimap.create();
        for (String line : lines) {
            if (!line.startsWith(HEADER_PREFIX)) {
                final TumorBAF region = fromString(line);
                result.put(region.chromosome(), region);
            }
        }
        return result;
    }

    @NotNull
    private static TumorBAF fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableTumorBAF.builder().chromosome(values[0]).position(Long.valueOf(values[1])).baf(Double.valueOf(values[2])).build();
    }
}

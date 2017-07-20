package com.hartwig.hmftools.common.purple.baf;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public enum TumorBAFFile {
    ;
    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "#";
    private static final String EXTENSION = ".purple.baf";

    @NotNull
    public static List<TumorBAF> read(@NotNull final String basePath, @NotNull final String sample) throws IOException {
        final String filePath = basePath + File.separator + sample + EXTENSION;
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String basePath, @NotNull final String sample, @NotNull List<TumorBAF> bafs)
            throws IOException {
        final String filePath = basePath + File.separator + sample + EXTENSION;
        Files.write(new File(filePath).toPath(), toLines(bafs));
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
    private static List<TumorBAF> fromLines(@NotNull List<String> lines) {
        return lines.stream().filter(x -> !x.startsWith(HEADER_PREFIX)).map(TumorBAFFile::fromString).collect(toList());
    }

    @NotNull
    private static TumorBAF fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableTumorBAF.builder().chromosome(values[0]).position(Long.valueOf(values[1])).baf(Double.valueOf(values[2])).build();
    }
}

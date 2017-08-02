package com.hartwig.hmftools.common.chromosome;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class ChromosomeLengthFile {

    private static final String DELIMITER = "\t";

    @NotNull
    public static List<ChromosomeLength> read(@NotNull final String filename) throws IOException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<ChromosomeLength> lengths) throws IOException {
        Files.write(new File(filename).toPath(), toLines(lengths));
    }

    @NotNull
    private static List<ChromosomeLength> fromLines(@NotNull List<String> lines) {
        return lines.stream().map(ChromosomeLengthFile::fromString).collect(Collectors.toList());
    }

    @NotNull
    static List<String> toLines(@NotNull final List<ChromosomeLength> lengths) {
        final List<String> lines = Lists.newArrayList();
        lengths.stream().map(ChromosomeLengthFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static ChromosomeLength fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableChromosomeLength.builder().chromosome(values[0]).position(Long.valueOf(values[2])).build();
    }

    @NotNull
    private static String toString(@NotNull final ChromosomeLength ratio) {
        return new StringJoiner(DELIMITER).add(String.valueOf(ratio.chromosome())).add(String.valueOf(ratio.position())).toString();
    }
}

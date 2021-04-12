package com.hartwig.hmftools.common.genome.chromosome;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class ChromosomeLengthFile {

    private static final String DELIMITER = "\t";
    private static final String EXTENSION = ".chr.len";

    private ChromosomeLengthFile() {
    }

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

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
    private static List<String> toLines(@NotNull final List<ChromosomeLength> lengths) {
        final List<String> lines = Lists.newArrayList();
        lengths.stream().map(ChromosomeLengthFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static ChromosomeLength fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableChromosomeLength.builder().chromosome(values[0]).length(Long.parseLong(values[1])).build();
    }

    @NotNull
    private static String toString(@NotNull final ChromosomeLength ratio) {
        return new StringJoiner(DELIMITER).add(ratio.chromosome()).add(String.valueOf(ratio.length())).toString();
    }
}

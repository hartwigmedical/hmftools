package com.hartwig.hmftools.common.cobalt;

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

@Deprecated
public enum ReadRatioFile {
    ;

    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "Chr";
    private static final String EXTENSION = ".cobalt.ratio";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static Multimap<String, ReadRatio> read(@NotNull final String filename) throws IOException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(@NotNull final String fileName, @NotNull Multimap<String, ReadRatio> ratios)
            throws IOException {
        List<ReadRatio> sorted = Lists.newArrayList(ratios.values());
        Collections.sort(sorted);
        write(fileName, sorted);
    }

    public static void write(@NotNull final String fileName, @NotNull List<ReadRatio> ratios)
            throws IOException {
        Files.write(new File(fileName).toPath(), toLines(ratios));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<ReadRatio> ratio) {
        final List<String> lines = Lists.newArrayList();
        lines.addAll(header());
        ratio.stream().map(ReadRatioFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static List<String> header() {
        final List<String> result = Lists.newArrayList();
        result.add(new StringJoiner(DELIMITER, "", "").add("Chromosome").add("Position").add("Ratio").toString());
        return result;
    }

    @NotNull
    private static String toString(@NotNull final ReadRatio ratio) {
        return new StringJoiner(DELIMITER).add(String.valueOf(ratio.chromosome()))
                .add(String.valueOf(ratio.position()))
                .add(String.valueOf(ratio.ratio()))
                .toString();
    }

    @NotNull
    private static Multimap<String, ReadRatio> fromLines(@NotNull List<String> lines) {
        Multimap<String, ReadRatio> result = ArrayListMultimap.create();
        for (String line : lines) {
            if (!line.startsWith(HEADER_PREFIX)) {
                final ReadRatio region = fromString(line);
                result.put(region.chromosome(), region);
            }
        }
        return result;
    }

    @NotNull
    private static ReadRatio fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableReadRatio.builder()
                .chromosome(values[0])
                .position(Long.valueOf(values[1]))
                .ratio(Double.valueOf(values[2]))
                .build();
    }
}

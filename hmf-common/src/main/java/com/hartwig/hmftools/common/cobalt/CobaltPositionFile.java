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

public class CobaltPositionFile {

    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "#";
    private static final String EXTENSION = ".cobalt";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(@NotNull final String fileName, @NotNull Multimap<String, CobaltPosition> ratios) throws IOException {
        List<CobaltPosition> sorted = Lists.newArrayList(ratios.values());
        Collections.sort(sorted);
        write(fileName, sorted);
    }

    public static void write(@NotNull final String fileName, @NotNull List<CobaltPosition> ratios) throws IOException {
        Files.write(new File(fileName).toPath(), toLines(ratios));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<CobaltPosition> ratio) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        ratio.stream().map(CobaltPositionFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, HEADER_PREFIX, "").add("Chromosome")
                .add("Position")
                .add("ReferenceReadCount")
                .add("ReferenceGCRatio")
                .add("ReferenceGCDiploidRatio")
                .add("TumorReadCount")
                .add("TumorGCRatio")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final CobaltPosition position) {
        return new StringJoiner(DELIMITER).add(String.valueOf(position.chromosome()))
                .add(String.valueOf(position.position()))
                .add(String.valueOf(position.referenceReadCount()))
                .add(String.valueOf(position.referenceGCRatio()))
                .add(String.valueOf(position.referenceGCDiploidRatio()))
                .add(String.valueOf(position.tumorReadCount()))
                .add(String.valueOf(position.tumorGCRatio()))
                .toString();
    }

    @NotNull
    private static Multimap<String, CobaltPosition> fromLines(@NotNull final List<String> lines) {

        final Multimap<String, CobaltPosition> result = ArrayListMultimap.create();
        for (String line : lines) {
            if (!line.startsWith(HEADER_PREFIX)) {
                final CobaltPosition position = fromLine(line);
                result.put(position.chromosome(), position);
            }
        }

        return result;
    }

    @NotNull
    private static CobaltPosition fromLine(@NotNull final String ratioLine) {
        final String[] values = ratioLine.split(DELIMITER);

        final String chromosome = values[0].trim();
        final long position = Long.valueOf(values[1].trim());

        return ImmutableCobaltPosition.builder()
                .chromosome(chromosome)
                .position(position)
                .referenceReadCount(Integer.valueOf(values[2].trim()))
                .referenceGCRatio(Double.valueOf(values[3].trim()))
                .referenceGCDiploidRatio(Double.valueOf(values[4].trim()))
                .tumorReadCount(Integer.valueOf(values[5].trim()))
                .tumorGCRatio(Double.valueOf(values[6].trim()))
                .build();
    }

}

package com.hartwig.hmftools.common.cobalt;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;

import org.jetbrains.annotations.NotNull;

public class CobaltRatioFile {

    private static final String DELIMITER = "\t";
    private static final String EXTENSION = ".cobalt";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static ListMultimap<String, CobaltRatio> read(@NotNull final String filename) throws IOException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(@NotNull final String fileName, @NotNull Multimap<Chromosome, CobaltRatio> ratios) throws IOException {
        List<CobaltRatio> sorted = Lists.newArrayList(ratios.values());
        Collections.sort(sorted);
        write(fileName, sorted);
    }

    private static void write(@NotNull final String fileName, @NotNull List<CobaltRatio> ratios) throws IOException {
        Files.write(new File(fileName).toPath(), toLines(ratios));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<CobaltRatio> ratio) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        ratio.stream().map(CobaltRatioFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("Chromosome")
                .add("Position")
                .add("ReferenceReadCount")
                .add("TumorReadCount")
                .add("ReferenceGCRatio")
                .add("TumorGCRatio")
                .add("ReferenceGCDiploidRatio")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final CobaltRatio position) {
        return new StringJoiner(DELIMITER).add(String.valueOf(position.chromosome()))
                .add(String.valueOf(position.position()))
                .add(String.valueOf(position.referenceReadCount()))
                .add(String.valueOf(position.tumorReadCount()))
                .add(String.valueOf(position.referenceGCRatio()))
                .add(String.valueOf(position.tumorGCRatio()))
                .add(String.valueOf(position.referenceGCDiploidRatio()))
                .toString();
    }

    @NotNull
    private static ListMultimap<String, CobaltRatio> fromLines(@NotNull final List<String> lines) {

        final ListMultimap<String, CobaltRatio> result = ArrayListMultimap.create();
        for (String line : lines) {
            if (!line.startsWith("Ch")) {
                final CobaltRatio position = fromLine(line);
                result.put(position.chromosome(), position);
            }
        }

        return result;
    }

    @NotNull
    static CobaltRatio fromLine(@NotNull final String ratioLine) {
        final String[] values = ratioLine.split(DELIMITER);

        final String chromosome = values[0].trim();
        final long position = Long.valueOf(values[1].trim());

        return ImmutableCobaltRatio.builder()
                .chromosome(chromosome)
                .position(position)
                .referenceReadCount(Integer.valueOf(values[2].trim()))
                .tumorReadCount(Integer.valueOf(values[3].trim()))
                .referenceGCRatio(Double.valueOf(values[4].trim()))
                .tumorGCRatio(Double.valueOf(values[5].trim()))
                .referenceGCDiploidRatio(Double.valueOf(values[6].trim()))
                .build();
    }

}

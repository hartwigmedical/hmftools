package com.hartwig.hmftools.common.cobalt;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class MedianRatioFile {

    private static final DecimalFormat FORMAT = new DecimalFormat("#.####");

    private static final String DELIMITER = "\t";
    private static final String EXTENSION = ".cobalt.ratio.median.tsv";

    private MedianRatioFile() {
    }

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static List<MedianRatio> read(@NotNull final String filename) throws IOException {
        return Files.readAllLines(new File(filename).toPath()).stream().skip(1).map(MedianRatioFile::fromLine).collect(Collectors.toList());
    }

    public static void write(@NotNull final String fileName, @NotNull List<MedianRatio> ratios) throws IOException {
        Files.write(new File(fileName).toPath(), toLines(ratios));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<MedianRatio> ratio) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        ratio.stream().map(MedianRatioFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("chromosome").add("medianRatio").add("count").toString();
    }

    @NotNull
    private static String toString(@NotNull final MedianRatio position) {
        return new StringJoiner(DELIMITER).add(position.chromosome())
                .add(String.valueOf(position.medianRatio()))
                .add(String.valueOf(position.count()))
                .toString();
    }

    @NotNull
    private static MedianRatio fromLine(@NotNull final String ratioLine) {
        final String[] values = ratioLine.split(DELIMITER);

        final String chromosome = values[0].trim();
        final double ratio = Double.parseDouble(values[1].trim());
        final int count = Integer.parseInt(values[2].trim());

        return ImmutableMedianRatio.builder().chromosome(chromosome).medianRatio(ratio).count(count).build();
    }
}

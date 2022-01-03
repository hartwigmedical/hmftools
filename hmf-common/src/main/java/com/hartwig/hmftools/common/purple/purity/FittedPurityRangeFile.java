package com.hartwig.hmftools.common.purple.purity;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.purple.purity.BestFit.bestFitPerPurity;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleCommon;

import org.jetbrains.annotations.NotNull;

public final class FittedPurityRangeFile {

    private static final DecimalFormat FORMAT = PurpleCommon.decimalFormat("0.0000");

    private static final String DELIMITER = "\t";
    private static final String COMMENT = "#";

    private static final String EXTENSION = ".purple.purity.range.tsv";
    private static final String EXTENSION_OLD = ".purple.purity.range";

    @NotNull
    public static String generateFilenameForReading(@NotNull final String basePath, @NotNull final String sample) {
        String filename = basePath + File.separator + sample + EXTENSION;
        return (new File(filename).exists()) ? filename : basePath + File.separator + sample + EXTENSION_OLD;
    }

    @NotNull
    public static List<FittedPurity> readAll(@NotNull final String basePath, @NotNull final String sample) throws IOException {
        final String filePath = generateFilenameForReading(basePath, sample);
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    @NotNull
    public static List<FittedPurity> readBestFitPerPurity(@NotNull final String basePath, @NotNull final String sample) throws IOException {
        return bestFitPerPurity(readAll(basePath, sample));
    }

    public static void write(@NotNull final String basePath, @NotNull final String sample, @NotNull final List<FittedPurity> purity)
            throws IOException {
        final String filePath = generateFilenameForWriting(basePath, sample);
        Files.write(new File(filePath).toPath(), toLines(purity));
    }

    @NotNull
    private static String generateFilenameForWriting(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull final List<FittedPurity> purity) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        purity.stream().map(FittedPurityRangeFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    @VisibleForTesting
    static List<FittedPurity> fromLines(@NotNull final List<String> lines) {
        return lines.stream()
                .filter(x -> !x.startsWith(COMMENT) && !x.startsWith("purity"))
                .map(FittedPurityRangeFile::fromString)
                .sorted()
                .collect(toList());
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "")
                .add("purity")
                .add("normFactor")
                .add("score")
                .add("diploidProportion")
                .add("ploidy")
                .add("somaticPenalty")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final FittedPurity purity) {
        return new StringJoiner(DELIMITER).add(FORMAT.format(purity.purity()))
                .add(FORMAT.format(purity.normFactor()))
                .add(FORMAT.format(purity.score()))
                .add(FORMAT.format(purity.diploidProportion()))
                .add(FORMAT.format(purity.ploidy()))
                .add(FORMAT.format(purity.somaticPenalty()))
                .toString();
    }

    @NotNull
    private static FittedPurity fromString(@NotNull final String purity) {
        String[] values = purity.split(DELIMITER);
        return ImmutableFittedPurity.builder()
                .purity(Double.parseDouble(values[0]))
                .normFactor(Double.parseDouble(values[1]))
                .score(Double.parseDouble(values[2]))
                .diploidProportion(Double.parseDouble(values[3]))
                .ploidy(Double.parseDouble(values[4]))
                .somaticPenalty(Double.parseDouble(values[5]))
                .build();
    }
}

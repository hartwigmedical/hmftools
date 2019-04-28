package com.hartwig.hmftools.common.purple.purity;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;
import java.util.TreeSet;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

public enum FittedPurityRangeFile {
    ;

    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final int MAX_RECORDS = 10000;
    private static final String DELIMITER = "\t";
    static final String HEADER_PREFIX = "#";
    private static final String EXTENSION = ".purple.purity.range";

    @NotNull
    public static List<FittedPurity> read(@NotNull final String basePath, @NotNull final String sample) throws IOException {
        final String filePath = basePath + File.separator + sample + EXTENSION;
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String basePath, @NotNull final String sample, @NotNull final List<FittedPurity> purity)
            throws IOException {
        final String filePath = basePath + File.separator + sample + EXTENSION;
        Files.write(new File(filePath).toPath(), toLines(purity));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<FittedPurity> purity) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        purity.stream().limit(MAX_RECORDS).map(FittedPurityRangeFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    static List<FittedPurity> fromLines(@NotNull final List<String> lines) {
        final List<FittedPurity> all =
                lines.stream().filter(x -> !x.startsWith(HEADER_PREFIX)).map(FittedPurityRangeFile::fromString).sorted().collect(toList());

        return bestFitPerPurity(all);
    }

    @NotNull
    static List<FittedPurity> bestFitPerPurity(@NotNull final List<FittedPurity> all) {
        Collections.sort(all);

        final List<FittedPurity> result = Lists.newArrayList();
        final TreeSet<Double> purities = new TreeSet<>(Doubles.comparator());
        for (FittedPurity fittedPurity : all) {
            if (purities.add(fittedPurity.purity())) {
                result.add(fittedPurity);
            }
        }

        return result;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, HEADER_PREFIX, "").add("Purity")
                .add("NormFactor")
                .add("Score")
                .add("DiploidProportion")
                .add("Ploidy")
                .add("SomaticPenalty")
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
                .purity(Double.valueOf(values[0]))
                .normFactor(Double.valueOf(values[1]))
                .score(Double.valueOf(values[2]))
                .diploidProportion(Double.valueOf(values[3]))
                .ploidy(Double.valueOf(values[4]))
                .somaticPenalty(Double.valueOf(values[5]))
                .build();
    }
}

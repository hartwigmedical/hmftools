package com.hartwig.hmftools.common.purple.purity;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public enum FittedPurityFile {
    ;

    private static final String DELIMITER = "\t";

    public static void writePurity(@NotNull final String filePath, @NotNull final List<FittedPurity> purity)
            throws IOException {
        Files.write(new File(filePath).toPath(), toLines(purity));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<FittedPurity> purity) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        purity.stream().limit(100).map(FittedPurityFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    static List<FittedPurity> fromLines(    @NotNull List<String> lines) {
        return lines.stream().filter(x -> !x.startsWith("#")).map(FittedPurityFile::fromString).collect(toList());
    }

    @NotNull
    private static String header() {
        return "#purity" + '\t' + "normFactor" + '\t' + "score" + '\t' + "modelBAFDeviation" + '\t'
                + "diploidProportion" + '\t' + "ploidy";
    }

    @NotNull
    private static String toString(@NotNull final FittedPurity purity) {
        return String.valueOf(purity.purity()) + DELIMITER + purity.normFactor() + DELIMITER + purity.score()
                + DELIMITER + purity.modelBAFDeviation() + DELIMITER + purity.diploidProportion() + DELIMITER
                + purity.ploidy();
    }

    @NotNull
    private static FittedPurity fromString(@NotNull final String purity) {
        String[] values = purity.split(DELIMITER);
        return ImmutableFittedPurity.builder()
                .purity(Double.valueOf(values[0]))
                .normFactor(Double.valueOf(values[1]))
                .score(Double.valueOf(values[2]))
                .modelBAFDeviation(Double.valueOf(values[3]))
                .diploidProportion(Double.valueOf(values[4]))
                .build();
    }
}

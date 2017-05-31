package com.hartwig.hmftools.common.purple.purity;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public enum FittedPurityScoreFile {
    ;

    private static final String DELIMITER = "\t";
    static final String HEADER_PREFIX = "#";

    @NotNull
    public static FittedPurityScore read(@NotNull final String filePath) throws IOException {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filePath, @NotNull final FittedPurityScore score)
            throws IOException {
        Files.write(new File(filePath).toPath(), toLines(score));
    }

    @NotNull
    static List<String> toLines(@NotNull final FittedPurityScore score) {
        return Lists.newArrayList(header(), toString(score));
    }

    @NotNull
    static FittedPurityScore fromLines(@NotNull List<String> lines) {
        assert (lines.size() > 1);
        return fromString(lines.get(1));
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, HEADER_PREFIX, "")
                .add("polyclonalProportion")
                .add("minPurity")
                .add("maxPurity")
                .add("minPloidy")
                .add("maxPloidy")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final FittedPurityScore score) {
        return new StringJoiner(DELIMITER).add(String.valueOf(score.polyclonalProportion()))
                .add(String.valueOf(score.minPurity()))
                .add(String.valueOf(score.maxPurity()))
                .add(String.valueOf(score.minPloidy()))
                .add(String.valueOf(score.maxPloidy()))
                .toString();
    }

    @NotNull
    private static FittedPurityScore fromString(@NotNull final String score) {
        String[] values = score.split(DELIMITER);
        return ImmutableFittedPurityScore.builder()
                .polyclonalProportion(Double.valueOf(values[0]))
                .minPurity(Double.valueOf(values[1]))
                .maxPurity(Double.valueOf(values[2]))
                .minPloidy(Double.valueOf(values[3]))
                .maxPloidy(Double.valueOf(values[4]))
                .build();
    }
}

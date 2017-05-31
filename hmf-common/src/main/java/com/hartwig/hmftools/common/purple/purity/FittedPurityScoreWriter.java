package com.hartwig.hmftools.common.purple.purity;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public enum FittedPurityScoreWriter {
    ;

    public static void writeScore(@NotNull final String filePath, @NotNull final FittedPurityScore score)
            throws IOException {
        final Collection<String> lines = Lists.newArrayList(header(), transform(score));
        Files.write(new File(filePath).toPath(), lines);
    }

    @NotNull
    private static String header() {
        return "polyclonalProportion" + '\t' + "minPurity" + '\t' + "maxPurity" + '\t' + "minPloidy" + '\t'
                + "maxPloidt";
    }

    @NotNull
    private static String transform(@NotNull final FittedPurityScore score) {
        return String.valueOf(score.polyclonalProportion()) + '\t' + score.minPurity() + '\t' + score.maxPurity()
                + '\t' + score.minPloidy() + '\t' + score.maxPloidy();
    }
}

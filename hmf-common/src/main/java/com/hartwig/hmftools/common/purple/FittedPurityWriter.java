package com.hartwig.hmftools.common.purple;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public enum FittedPurityWriter {
    ;

    public static void writePurity(@NotNull final String filePath, @NotNull final Collection<FittedPurity> purity)
            throws IOException {
        final Collection<String> lines = Lists.newArrayList();
        lines.add(header());
        purity.stream().limit(100).map(FittedPurityWriter::transform).forEach(lines::add);

        Files.write(new File(filePath).toPath(), lines);
    }

    @NotNull
    private static String header() {
        return "purity" + '\t' + "normFactor" + '\t' + "score" + '\t' + "modelBAFDeviation" + '\t'
                + "diploidProportion" + '\t' + "ploidy";
    }

    @NotNull
    private static String transform(@NotNull final FittedPurity purity) {
        return String.valueOf(purity.purity()) + '\t' + purity.normFactor() + '\t' + purity.score() + '\t'
                + purity.modelBAFDeviation() + '\t' + purity.diploidProportion() + '\t' + purity.ploidy();
    }
}

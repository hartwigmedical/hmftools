package com.hartwig.hmftools.common.purple;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public enum FittedPurityWriter {
    ;

    public static void writePurity(@NotNull final String filePath, @NotNull Collection<FittedPurity> purity) throws IOException {

        final Collection<String> lines = Lists.newArrayList();
        lines.add(header());
        purity.stream().limit(100).map(FittedPurityWriter::transform).forEach(lines::add);

        Files.write(new File(filePath).toPath(), lines);
    }

    private static String header() {
        return new StringBuilder()
                .append("purity").append('\t')
                .append("normFactor").append('\t')
                .append("score").append('\t')
                .append("modelBAFDeviation").append('\t')
                .append("diplodProportion").append('\t')
                .append("ploidy")
                .toString();
    }

    private static String transform(FittedPurity purity) {
        return new StringBuilder()
                .append(purity.purity()).append('\t')
                .append(purity.normFactor()).append('\t')
                .append(purity.score()).append('\t')
                .append(purity.modelBAFDeviation()).append('\t')
                .append(purity.diplodProportion()).append('\t')
                .append(purity.ploidy())
                .toString();
    }

}

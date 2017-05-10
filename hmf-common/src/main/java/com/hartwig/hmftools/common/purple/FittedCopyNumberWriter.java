package com.hartwig.hmftools.common.purple;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;

public enum FittedCopyNumberWriter {
    ;

    public static void writeCopyNumber(@NotNull final String filePath, @NotNull Collection<FittedCopyNumber> copyNumbers) throws IOException {

        final Collection<String> lines = Lists.newArrayList();
        lines.add(header());
        copyNumbers.stream().map(FittedCopyNumberWriter::transform).forEach(lines::add);

        Files.write(new File(filePath).toPath(), lines);
    }

    private static String header() {
        return new StringBuilder()
                .append("chromosome").append('\t')
                .append("start").append('\t')
                .append("end").append('\t')
                .append("copyNumber").append('\t')
                .append("alteration").append('\t')
                .append("fittedPloidy").append('\t')
                .append("actualBAF").append('\t')
                .append("modelBAF").append('\t')
                .append("bafDeviation").append('\t')
                .append("actualCNVRatio").append('\t')
                .append("modelCNVRatio").append('\t')
                .append("cnvDeviation").append('\t')
                .append("deviation")
                .toString();
    }

    private static String transform(FittedCopyNumber copyNumber) {
        return new StringBuilder()
                .append(copyNumber.chromosome()).append('\t')
                .append(copyNumber.start()).append('\t')
                .append(copyNumber.end()).append('\t')
                .append(copyNumber.value()).append('\t')
                .append(alteration(copyNumber)).append('\t')
                .append(copyNumber.fittedPloidy()).append('\t')
                .append(copyNumber.actualBAF()).append('\t')
                .append(copyNumber.modelBAF()).append('\t')
                .append(copyNumber.bafDeviation()).append('\t')
                .append(copyNumber.actualCNVRatio()).append('\t')
                .append(copyNumber.modelCNVRatio()).append('\t')
                .append(copyNumber.cnvDeviation()).append('\t')
                .append(copyNumber.deviation())
                .toString();

    }

    private static String alteration(CopyNumber copyNumber) {
        return copyNumber.isGain() ? "gain" : copyNumber.isLoss() ? "lost" : "neutral";
    }
}

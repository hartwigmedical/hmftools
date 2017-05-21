package com.hartwig.hmftools.common.purple;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.CopyNumber;

import org.jetbrains.annotations.NotNull;

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
                .append("status").append('\t')
                .append("fittedPloidy").append('\t')
                .append("bafCount").append('\t')
                .append("observedBAF").append('\t')
                .append("modelBAF").append('\t')
                .append("bafDeviation").append('\t')
                .append("observedTumorRatio").append('\t')
                .append("normalRatio").append('\t')
                .append("modelTumorRatio").append('\t')
                .append("cnvDeviation").append('\t')
                .append("deviation").append('\t')
                .append("ratioOfRatios").append('\t')
                .append("broadRatioOfRatios").append('\t')
                .append("broadBAF").append('\t')
                .append("segmentRatio").append('\t')
                .append("segmentBAF")
                .toString();
    }

    private static String transform(FittedCopyNumber copyNumber) {
        return new StringBuilder()
                .append(copyNumber.chromosome()).append('\t')
                .append(copyNumber.start()).append('\t')
                .append(copyNumber.end()).append('\t')
                .append(copyNumber.value()).append('\t')
                .append(alteration(copyNumber)).append('\t')
                .append(copyNumber.status().toString().toLowerCase()).append('\t')
                .append(copyNumber.fittedPloidy()).append('\t')
                .append(copyNumber.bafCount()).append('\t')
                .append(copyNumber.observedBAF()).append('\t')
                .append(copyNumber.modelBAF()).append('\t')
                .append(copyNumber.bafDeviation()).append('\t')
                .append(copyNumber.observedTumorRatio()).append('\t')
                .append(copyNumber.normalRatio()).append('\t')
                .append(copyNumber.modelTumorRatio()).append('\t')
                .append(copyNumber.cnvDeviation()).append('\t')
                .append(copyNumber.deviation()).append('\t')
                .append(copyNumber.ratioOfRatios()).append('\t')
                .append(copyNumber.broadRatioOfRatios()).append('\t')
                .append(copyNumber.broadBAF()).append('\t')
                .append(copyNumber.segmentRatioOfRatios()).append('\t')
                .append(copyNumber.segmentBAF())
                .toString();

    }

    private static String alteration(CopyNumber copyNumber) {
        return copyNumber.isGain() ? "gain" : copyNumber.isLoss() ? "lost" : "neutral";
    }
}

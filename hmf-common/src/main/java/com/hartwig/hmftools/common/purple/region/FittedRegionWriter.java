package com.hartwig.hmftools.common.purple.region;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.CopyNumber;

import org.jetbrains.annotations.NotNull;

public enum FittedRegionWriter {
    ;

    private static final String EXTENSION = ".purple.fitted";


    public static void writeCopyNumber(@NotNull final String basePath, @NotNull final String sample,
            @NotNull Collection<FittedRegion> copyNumbers) throws IOException {
        final String filePath = basePath + File.separator + sample + EXTENSION;
        final Collection<String> lines = Lists.newArrayList();
        lines.add(header());
        copyNumbers.stream().map(FittedRegionWriter::transform).forEach(lines::add);

        Files.write(new File(filePath).toPath(), lines);
    }

    @NotNull
    private static String header() {
        return "chromosome" + '\t' + "start" + '\t' + "end" + '\t'
                + "status" + '\t' + "fittedPloidy" + '\t' + "bafCount" + '\t' + "observedBAF" + '\t'
                + "purityAdjustedBAF" + '\t' + "modelBAF" + '\t' + "bafDeviation" + '\t' + "observedTumorRatio" + '\t'
                + "observedNormalRatio" + '\t' + "modelTumorRatio" + '\t' + "cnvDeviation" + '\t' + "deviation" + '\t'
                + "tumorCopyNumber" + '\t' + "broadTumorCopyNumber" + '\t' + "broadBAF" + '\t'
                + "segmentTumorCopyNumber" + '\t' + "segmentBAF";
    }

    @NotNull
    private static String transform(@NotNull final FittedRegion copyNumber) {
        return copyNumber.chromosome() + '\t' + copyNumber.start() + '\t' + copyNumber.end() + '\t'
                + copyNumber.status().toString().toLowerCase() + '\t' + copyNumber.fittedPloidy() + '\t'
                + copyNumber.bafCount() + '\t' + copyNumber.observedBAF() + '\t' + copyNumber.purityAdjustedBAF()
                + '\t' + copyNumber.modelBAF() + '\t' + copyNumber.bafDeviation() + '\t'
                + copyNumber.observedTumorRatio() + '\t' + copyNumber.observedNormalRatio() + '\t'
                + copyNumber.modelTumorRatio() + '\t' + copyNumber.cnvDeviation() + '\t' + copyNumber.deviation()
                + '\t' + copyNumber.tumorCopyNumber() + '\t' + copyNumber.broadTumorCopyNumber() + '\t'
                + copyNumber.broadBAF() + '\t' + copyNumber.segmentTumorCopyNumber() + '\t' + copyNumber.segmentBAF();
    }

    private static String alteration(@NotNull final CopyNumber copyNumber) {
        return copyNumber.isGain() ? "gain" : copyNumber.isLoss() ? "lost" : "neutral";
    }
}

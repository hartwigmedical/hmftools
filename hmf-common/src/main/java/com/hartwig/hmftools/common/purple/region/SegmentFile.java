package com.hartwig.hmftools.common.purple.region;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class SegmentFile {
    
    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final String EXTENSION = ".purple.segment.tsv";
    private static final String DELIMITER = "\t";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(@NotNull final String filePath, @NotNull Collection<FittedRegion> fittedRegions) throws IOException {
        final Collection<String> lines = Lists.newArrayList();
        lines.add(header());
        fittedRegions.stream().map(SegmentFile::toString).forEach(lines::add);
        Files.write(new File(filePath).toPath(), lines);
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("chromosome")
                .add("start")
                .add("end")
                .add("germlineStatus")
                .add("bafCount")
                .add("observedBAF")
                .add("minorAlleleCopyNumber")
                .add("minorAlleleCopyNumberDeviation")
                .add("observedTumorRatio")
                .add("observedNormalRatio")
                .add("unnormalisedObservedNormalRatio")
                .add("majorAlleleCopyNumber")
                .add("majorAlleleCopyNumberDeviation")
                .add("deviationPenalty")
                .add("tumorCopyNumber")
                .add("fittedTumorCopyNumber")
                .add("fittedBAF")
                .add("refNormalisedCopyNumber")
                .add("ratioSupport")
                .add("support")
                .add("depthWindowCount")
                .add("tumorBAF")
                .add("gcContent")
                .add("svCluster")
                .add("eventPenalty")
                .add("minStart")
                .add("maxStart")
                .toString();
    }

    @NotNull
    static String toString(@NotNull final FittedRegion copyNumber) {
        return new StringJoiner(DELIMITER).add(copyNumber.chromosome())
                .add(String.valueOf(copyNumber.start()))
                .add(String.valueOf(copyNumber.end()))
                .add(String.valueOf(copyNumber.status()))
                .add(String.valueOf(copyNumber.bafCount()))
                .add(FORMAT.format(copyNumber.observedBAF()))
                .add(FORMAT.format(copyNumber.minorAlleleCopyNumber()))
                .add(FORMAT.format(copyNumber.minorAlleleCopyNumberDeviation()))
                .add(FORMAT.format(copyNumber.observedTumorRatio()))
                .add(FORMAT.format(copyNumber.observedNormalRatio()))
                .add(FORMAT.format(copyNumber.unnormalisedObservedNormalRatio()))
                .add(FORMAT.format(copyNumber.majorAlleleCopyNumber()))
                .add(FORMAT.format(copyNumber.majorAlleleCopyNumberDeviation()))
                .add(FORMAT.format(copyNumber.deviationPenalty()))
                .add(FORMAT.format(copyNumber.tumorCopyNumber()))
                .add(FORMAT.format(copyNumber.fittedTumorCopyNumber()))
                .add(FORMAT.format(copyNumber.fittedBAF()))
                .add(FORMAT.format(copyNumber.refNormalisedCopyNumber()))
                .add(String.valueOf(copyNumber.ratioSupport()))
                .add(String.valueOf(copyNumber.support()))
                .add(String.valueOf(copyNumber.depthWindowCount()))
                .add(FORMAT.format(copyNumber.tumorBAF()))
                .add(FORMAT.format(copyNumber.gcContent()))
                .add(String.valueOf(copyNumber.svCluster()))
                .add(FORMAT.format(copyNumber.eventPenalty()))
                .add(String.valueOf(copyNumber.minStart()))
                .add(String.valueOf(copyNumber.maxStart()))
                .toString();
    }
}

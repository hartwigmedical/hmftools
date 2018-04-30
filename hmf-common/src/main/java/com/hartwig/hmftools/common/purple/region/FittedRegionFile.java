package com.hartwig.hmftools.common.purple.region;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

public enum FittedRegionFile {
    ;

    private static final String EXTENSION = ".purple.fitted";
    private static final String DELIMITER = "\t";
    static final String HEADER_PREFIX = "#";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(@NotNull final String filePath, @NotNull Collection<FittedRegion> copyNumbers) throws IOException {
        final Collection<String> lines = Lists.newArrayList();
        lines.add(header());
        copyNumbers.stream().map(FittedRegionFile::toString).forEach(lines::add);
        Files.write(new File(filePath).toPath(), lines);
    }

    @NotNull
    public static List<FittedRegion> read(final String filePath) throws IOException {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    @NotNull
    static List<FittedRegion> fromLines(@NotNull List<String> lines) {
        return lines.stream().filter(x -> !x.startsWith(HEADER_PREFIX)).map(FittedRegionFile::fromString).collect(toList());
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, HEADER_PREFIX, "").add("chromosome")
                .add("start")
                .add("end")
                .add("germlineStatus")
                .add("fittedPloidy")
                .add("bafCount")
                .add("observedBAF")
                .add("modelBAF")
                .add("bafDeviation")
                .add("observedTumorRatio")
                .add("observedNormalRatio")
                .add("modelTumorRatio")
                .add("cnvDeviation")
                .add("deviation")
                .add("tumorCopyNumber")
                .add("segmentTumorCopyNumber")
                .add("segmentBAF")
                .add("refNormalisedCopyNumber")
                .add("ratioSupport")
                .add("support")
                .add("observedTumorRatioCount")
                .add("tumorBAF")
                .add("gcContent")
                .add("svCluster")
                .add("ploidyPenalty")
                .toString();
    }

    @NotNull
    static String toString(@NotNull final FittedRegion copyNumber) {
        return new StringJoiner(DELIMITER).add(String.valueOf(copyNumber.chromosome()))
                .add(String.valueOf(copyNumber.start()))
                .add(String.valueOf(copyNumber.end()))
                .add(String.valueOf(copyNumber.status()))
                .add(String.valueOf(copyNumber.modelPloidy()))
                .add(String.valueOf(copyNumber.bafCount()))
                .add(String.valueOf(copyNumber.observedBAF()))
                .add(String.valueOf(copyNumber.modelBAF()))
                .add(String.valueOf(copyNumber.bafDeviation()))
                .add(String.valueOf(copyNumber.observedTumorRatio()))
                .add(String.valueOf(copyNumber.observedNormalRatio()))
                .add(String.valueOf(copyNumber.modelTumorRatio()))
                .add(String.valueOf(copyNumber.cnvDeviation()))
                .add(String.valueOf(copyNumber.deviation()))
                .add(String.valueOf(copyNumber.tumorCopyNumber()))
                .add(String.valueOf(copyNumber.segmentTumorCopyNumber()))
                .add(String.valueOf(copyNumber.segmentBAF()))
                .add(String.valueOf(copyNumber.refNormalisedCopyNumber()))
                .add(String.valueOf(copyNumber.ratioSupport()))
                .add(String.valueOf(copyNumber.support()))
                .add(String.valueOf(copyNumber.observedTumorRatioCount()))
                .add(String.valueOf(copyNumber.tumorBAF()))
                .add(String.valueOf(copyNumber.gcContent()))
                .add(String.valueOf(copyNumber.svCluster()))
                .add(String.valueOf(copyNumber.ploidyPenalty()))
                .toString();
    }

    @NotNull
    static FittedRegion fromString(@NotNull final String region) {
        String[] values = region.split(DELIMITER);
        ImmutableFittedRegion.Builder builder = ImmutableFittedRegion.builder()
                .chromosome(values[0])
                .start(Long.valueOf(values[1]))
                .end(Long.valueOf(values[2]))
                .status(GermlineStatus.fromString(values[3]))
                .modelPloidy(Integer.valueOf(values[4]))
                .bafCount(Integer.valueOf(values[5]))
                .observedBAF(Double.valueOf(values[6]))
                .modelBAF(Double.valueOf(values[7]))
                .bafDeviation(Double.valueOf(values[8]))
                .observedTumorRatio(Double.valueOf(values[9]))
                .observedNormalRatio(Double.valueOf(values[10]))
                .modelTumorRatio(Double.valueOf(values[11]))
                .cnvDeviation(Double.valueOf(values[12]))
                .deviation(Double.valueOf(values[13]))
                .tumorCopyNumber(Double.valueOf(values[14]))
                .segmentTumorCopyNumber(Double.valueOf(values[15]))
                .segmentBAF(Double.valueOf(values[16]))
                .refNormalisedCopyNumber(Double.valueOf(values[17]))
                .ratioSupport(Boolean.valueOf(values[18]))
                .support(SegmentSupport.valueOf(values[19]))
                .observedTumorRatioCount(Integer.valueOf(values[20]))
                .tumorBAF(Double.valueOf(values[21]))
                .gcContent(Double.valueOf(values[22]))
                .svCluster(Boolean.valueOf(values[23]))
                .ploidyPenalty(0);

        if (values.length == 25) {
            builder.ploidyPenalty(Double.valueOf(values[24]));
        }

        return builder.build();
    }
}

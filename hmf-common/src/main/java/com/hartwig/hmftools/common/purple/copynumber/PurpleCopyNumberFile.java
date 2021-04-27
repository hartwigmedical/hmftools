package com.hartwig.hmftools.common.purple.copynumber;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

public final class PurpleCopyNumberFile {

    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final String DELIMITER = "\t";

    private static final String SOMATIC_EXTENSION = ".purple.cnv.somatic.tsv";
    private static final String GERMLINE_EXTENSION = ".purple.cnv.germline.tsv";

    private static final String SOMATIC_EXTENSION_OLD = ".purple.cnv";
    private static final String GERMLINE_EXTENSION_OLD = ".purple.germline.cnv";

    private PurpleCopyNumberFile() {
    }

    @NotNull
    public static String generateFilenameForWriting(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + SOMATIC_EXTENSION;
    }

    @NotNull
    public static String generateFilenameForReading(@NotNull final String basePath, @NotNull final String sample) {
        String filename = basePath + File.separator + sample + SOMATIC_EXTENSION;
        return (new File(filename).exists()) ? filename : basePath + File.separator + sample + SOMATIC_EXTENSION_OLD;
    }

    @NotNull
    public static String generateGermlineFilenameForWriting(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + GERMLINE_EXTENSION;
    }

    @NotNull
    public static String generateGermlineFilenameForReading(@NotNull final String basePath, @NotNull final String sample) {
        String filename = basePath + File.separator + sample + GERMLINE_EXTENSION;
        return (new File(filename).exists()) ? filename : basePath + File.separator + sample + GERMLINE_EXTENSION_OLD;
    }

    @NotNull
    public static List<PurpleCopyNumber> read(final String filePath) throws IOException {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<PurpleCopyNumber> copyNumbers) throws IOException {
        Files.write(new File(filename).toPath(), toLines(copyNumbers));
    }

    @VisibleForTesting
    @NotNull
    public static List<String> toLines(@NotNull final List<PurpleCopyNumber> copyNumbers) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        copyNumbers.stream().map(PurpleCopyNumberFile::toString).forEach(lines::add);
        return lines;
    }

    @VisibleForTesting
    @NotNull
    public static List<PurpleCopyNumber> fromLines(@NotNull List<String> lines) {
        return lines.stream()
                .skip(1)
                .map(PurpleCopyNumberFile::fromString)
                .collect(toList());
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("chromosome")
                .add("start")
                .add("end")
                .add("copyNumber")
                .add("bafCount")
                .add("observedBAF")
                .add("baf")
                .add("segmentStartSupport")
                .add("segmentEndSupport")
                .add("method")
                .add("depthWindowCount")
                .add("gcContent")
                .add("minStart")
                .add("maxStart")
                .add("minorAlleleCopyNumber")
                .add("majorAlleleCopyNumber")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final PurpleCopyNumber copyNumber) {
        return new StringJoiner(DELIMITER).add(copyNumber.chromosome())
                .add(String.valueOf(copyNumber.start()))
                .add(String.valueOf(copyNumber.end()))
                .add(FORMAT.format(copyNumber.averageTumorCopyNumber()))
                .add(String.valueOf(copyNumber.bafCount()))
                .add(FORMAT.format(copyNumber.averageObservedBAF()))
                .add(FORMAT.format(copyNumber.averageActualBAF()))
                .add(String.valueOf(copyNumber.segmentStartSupport()))
                .add(String.valueOf(copyNumber.segmentEndSupport()))
                .add(String.valueOf(copyNumber.method()))
                .add(String.valueOf(copyNumber.depthWindowCount()))
                .add(FORMAT.format(copyNumber.gcContent()))
                .add(String.valueOf(copyNumber.minStart()))
                .add(String.valueOf(copyNumber.maxStart()))
                .add(FORMAT.format(copyNumber.minorAlleleCopyNumber()))
                .add(FORMAT.format(copyNumber.majorAlleleCopyNumber()))
                .toString();
    }

    @NotNull
    private static PurpleCopyNumber fromString(@NotNull final String copyNumber) {
        String[] values = copyNumber.split(DELIMITER);
        final ImmutablePurpleCopyNumber.Builder builder = ImmutablePurpleCopyNumber.builder()
                .chromosome(values[0])
                .start(Long.parseLong(values[1]))
                .end(Long.parseLong(values[2]))
                .averageTumorCopyNumber(Double.parseDouble(values[3]))
                .bafCount(Integer.parseInt(values[4]))
                .averageObservedBAF(Double.parseDouble(values[5]))
                .averageActualBAF(Double.parseDouble(values[6]))
                .segmentStartSupport(SegmentSupport.valueOf(values[7]))
                .segmentEndSupport(SegmentSupport.valueOf(values[8]))
                .method(CopyNumberMethod.valueOf(values[9]))
                .depthWindowCount(Integer.parseInt(values[10]))
                .gcContent(Double.parseDouble(values[11]))
                .minStart(Long.parseLong(values[12]))
                .maxStart(Long.parseLong(values[13]));

        return builder.build();
    }
}

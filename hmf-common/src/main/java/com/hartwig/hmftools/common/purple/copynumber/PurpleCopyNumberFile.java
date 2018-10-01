package com.hartwig.hmftools.common.purple.copynumber;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

public enum PurpleCopyNumberFile {
    ;

    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final String DELIMITER = "\t";
    static final String HEADER_PREFIX = "#";
    private static final String EXTENSION = ".purple.cnv";
    private static final String GERMLINE_EXTENSION = ".purple.germline.cnv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static String generateGermlineFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + GERMLINE_EXTENSION;
    }

    @NotNull
    public static List<PurpleCopyNumber> read(final String filePath) throws IOException {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<PurpleCopyNumber> copyNumbers) throws IOException {
        Files.write(new File(filename).toPath(), toLines(copyNumbers));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<PurpleCopyNumber> purity) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        purity.stream().map(PurpleCopyNumberFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    static List<PurpleCopyNumber> fromLines(@NotNull List<String> lines) {
        return lines.stream().filter(x -> !x.startsWith(HEADER_PREFIX)).map(PurpleCopyNumberFile::fromString).collect(toList());
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, HEADER_PREFIX, "").add("chromosome")
                .add("start")
                .add("end")
                .add("copyNumber")
                .add("bafCount")
                .add("observedBAF")
                .add("actualBAF")
                .add("segmentStartSupport")
                .add("segmentEndSupport")
                .add("method")
                .add("depthWindowCount")
                .add("gcContent")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final PurpleCopyNumber copyNumber) {
        return new StringJoiner(DELIMITER).add(String.valueOf(copyNumber.chromosome()))
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
                .toString();
    }

    @NotNull
    private static PurpleCopyNumber fromString(@NotNull final String copyNumber) {
        String[] values = copyNumber.split(DELIMITER);
        final ImmutablePurpleCopyNumber.Builder builder = ImmutablePurpleCopyNumber.builder()
                .chromosome(values[0])
                .start(Long.valueOf(values[1]))
                .end(Long.valueOf(values[2]))
                .averageTumorCopyNumber(Double.valueOf(values[3]))
                .bafCount(Integer.valueOf(values[4]))
                .averageObservedBAF(Double.valueOf(values[5]))
                .averageActualBAF(Double.valueOf(values[6]))
                .segmentStartSupport(SegmentSupport.valueOf(values[7]))
                .segmentEndSupport(SegmentSupport.valueOf(values[8]))
                .method(CopyNumberMethod.valueOf(values[9]))
                .depthWindowCount(0)
                .gcContent(0);
        if (values.length > 10) {
            builder.depthWindowCount(Integer.valueOf(values[10]));
        }
        if (values.length > 11) {
            builder.gcContent(Double.valueOf(values[11]));
        }

        return builder.build();
    }
}

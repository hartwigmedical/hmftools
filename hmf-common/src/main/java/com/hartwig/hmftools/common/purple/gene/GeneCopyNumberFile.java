package com.hartwig.hmftools.common.purple.gene;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

public final class GeneCopyNumberFile {

    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final String DELIMITER = "\t";

    private static final String EXTENSION = ".purple.cnv.gene.tsv";
    private static final String EXTENSION_OLD = ".purple.gene.cnv";

    private GeneCopyNumberFile() {
    }

    @NotNull
    public static String generateFilenameForWriting(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static String generateFilenameForReading(@NotNull final String basePath, @NotNull final String sample) {
        String filename = basePath + File.separator + sample + EXTENSION;
        return (new File(filename).exists()) ? filename : basePath + File.separator + sample + EXTENSION_OLD;
    }

    @NotNull
    public static List<GeneCopyNumber> read(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    public static void write(@NotNull final String fileName, @NotNull List<GeneCopyNumber> geneCopyNumbers) throws IOException {
        Files.write(new File(fileName).toPath(), toLines(geneCopyNumbers));
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull final List<GeneCopyNumber> ratio) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        ratio.stream().map(GeneCopyNumberFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("chromosome")
                .add("start")
                .add("end")
                .add("gene")
                .add("minCopyNumber")
                .add("maxCopyNumber")
                .add("unused")
                .add("somaticRegions")
                .add("germlineHomDeletionRegions")
                .add("germlineHetToHomDeletionRegions")
                .add("transcriptId")
                .add("transcriptVersion")
                .add("chromosomeBand")
                .add("minRegions")
                .add("minRegionStart")
                .add("minRegionEnd")
                .add("minRegionStartSupport")
                .add("minRegionEndSupport")
                .add("minRegionMethod")
                .add("minMinorAlleleCopyNumber")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final GeneCopyNumber geneCopyNumber) {

        return new StringJoiner(DELIMITER).add(geneCopyNumber.chromosome())
                .add(String.valueOf(geneCopyNumber.start()))
                .add(String.valueOf(geneCopyNumber.end()))
                .add(geneCopyNumber.gene())
                .add(FORMAT.format(geneCopyNumber.minCopyNumber()))
                .add(FORMAT.format(geneCopyNumber.maxCopyNumber()))
                .add(String.valueOf(0))
                .add(String.valueOf(geneCopyNumber.somaticRegions()))
                .add(String.valueOf(geneCopyNumber.germlineHomRegions()))
                .add(String.valueOf(geneCopyNumber.germlineHet2HomRegions()))
                .add(geneCopyNumber.transcriptID())
                .add(String.valueOf(geneCopyNumber.transcriptVersion()))
                .add(geneCopyNumber.chromosomeBand())
                .add(String.valueOf(geneCopyNumber.minRegions()))
                .add(String.valueOf(geneCopyNumber.minRegionStart()))
                .add(String.valueOf(geneCopyNumber.minRegionEnd()))
                .add(String.valueOf(geneCopyNumber.minRegionStartSupport()))
                .add(String.valueOf(geneCopyNumber.minRegionEndSupport()))
                .add(String.valueOf(geneCopyNumber.minRegionMethod()))
                .add(FORMAT.format(geneCopyNumber.minMinorAlleleCopyNumber()))
                .toString();
    }

    @NotNull
    @VisibleForTesting
    static List<GeneCopyNumber> fromLines(@NotNull List<String> lines) {
        return lines.stream()
                .skip(1)
                .map(GeneCopyNumberFile::fromString)
                .collect(toList());
    }

    @VisibleForTesting
    @NotNull
    static GeneCopyNumber fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);

        final ImmutableGeneCopyNumber.Builder builder = ImmutableGeneCopyNumber.builder()
                .chromosome(values[0])
                .start(Long.parseLong(values[1]))
                .end(Long.parseLong(values[2]))
                .gene(values[3])
                .minCopyNumber(Double.parseDouble(values[4]))
                .maxCopyNumber(Double.parseDouble(values[5]))
                .somaticRegions(Integer.parseInt(values[7]))
                .germlineHomRegions(Integer.parseInt(values[8]))
                .germlineHet2HomRegions(Integer.parseInt(values[9]))
                .transcriptID(values[10])
                .transcriptVersion(Integer.parseInt(values[11]))
                .chromosomeBand(values[12]);

        builder.minRegions(0)
                .minRegionStart(Long.parseLong(values[1]))
                .minRegionEnd(Long.parseLong(values[2]))
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEndSupport(SegmentSupport.NONE)
                .minMinorAlleleCopyNumber(0);

        if (values.length >= 19) {
            builder.minRegions(Integer.parseInt(values[13]))
                    .minRegionStart(Long.parseLong(values[14]))
                    .minRegionEnd(Long.parseLong(values[15]))
                    .minRegionStartSupport(SegmentSupport.valueOf(values[16]))
                    .minRegionEndSupport(SegmentSupport.valueOf(values[17]))
                    .minRegionMethod(CopyNumberMethod.valueOf(values[18]));
        }

        builder.minMinorAlleleCopyNumber(Double.parseDouble(values[values.length - 1]));

        return builder.build();
    }
}

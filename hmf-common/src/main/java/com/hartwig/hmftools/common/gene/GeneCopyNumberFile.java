package com.hartwig.hmftools.common.gene;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

public class GeneCopyNumberFile {
    ;

    static final String DELIMITER = "\t";
    static final String HEADER_PREFIX = "Chr";
    private static final String EXTENSION = ".purple.gene.cnv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static List<GeneCopyNumber> read(@NotNull final String fileName) throws IOException {
        return read(new File(fileName));
    }

    @NotNull
    public static List<GeneCopyNumber> read(@NotNull final File file) throws IOException {
        return fromLines(Files.readAllLines(file.toPath()));
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
        return new StringJoiner(DELIMITER, "", "").add("Chromosome")
                .add("Start")
                .add("End")
                .add("Gene")
                .add("MinCopyNumber")
                .add("MaxCopyNumber")
                .add("MeanCopyNumber")
                .add("SomaticRegions")
                .add("GermlineHomRegions")
                .add("GermlineHet2HomRegions")
                .add("TranscriptId")
                .add("TranscriptVersion")
                .add("ChromosomeBand")
                .add("MinRegions")
                .add("MinRegionStart")
                .add("MinRegionEnd")
                .add("MinRegionStartSupport")
                .add("MinRegionEndSupport")
                .add("MinRegionMethod")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final GeneCopyNumber geneCopyNumber) {
        return new StringJoiner(DELIMITER).add(String.valueOf(geneCopyNumber.chromosome()))
                .add(String.valueOf(geneCopyNumber.start()))
                .add(String.valueOf(geneCopyNumber.end()))
                .add(String.valueOf(geneCopyNumber.gene()))
                .add(String.valueOf(geneCopyNumber.minCopyNumber()))
                .add(String.valueOf(geneCopyNumber.maxCopyNumber()))
                .add(String.valueOf(geneCopyNumber.meanCopyNumber()))
                .add(String.valueOf(geneCopyNumber.somaticRegions()))
                .add(String.valueOf(geneCopyNumber.germlineHomRegions()))
                .add(String.valueOf(geneCopyNumber.germlineHet2HomRegions()))
                .add(String.valueOf(geneCopyNumber.transcriptID()))
                .add(String.valueOf(geneCopyNumber.transcriptVersion()))
                .add(String.valueOf(geneCopyNumber.chromosomeBand()))
                .add(String.valueOf(geneCopyNumber.minRegions()))
                .add(String.valueOf(geneCopyNumber.minRegionStart()))
                .add(String.valueOf(geneCopyNumber.minRegionEnd()))
                .add(String.valueOf(geneCopyNumber.minRegionStartSupport()))
                .add(String.valueOf(geneCopyNumber.minRegionEndSupport()))
                .add(String.valueOf(geneCopyNumber.minRegionMethod()))
                .toString();
    }

    @NotNull
    @VisibleForTesting
    static List<GeneCopyNumber> fromLines(@NotNull List<String> lines) {
        return lines.stream().filter(x -> !x.startsWith(HEADER_PREFIX)).map(GeneCopyNumberFile::fromString).collect(toList());
    }

    @VisibleForTesting
    @NotNull
    static GeneCopyNumber fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);

        final ImmutableGeneCopyNumber.Builder builder = ImmutableGeneCopyNumber.builder()
                .chromosome(values[0])
                .start(Long.valueOf(values[1]))
                .end(Long.valueOf(values[2]))
                .gene(values[3])
                .minCopyNumber(Double.valueOf(values[4]))
                .maxCopyNumber(Double.valueOf(values[5]))
                .meanCopyNumber(Double.valueOf(values[6]))
                .somaticRegions(Integer.valueOf(values[7]))
                .germlineHomRegions(Integer.valueOf(values[8]))
                .germlineHet2HomRegions(Integer.valueOf(values[9]))
                .transcriptID(values[10])
                .transcriptVersion(Integer.valueOf(values[11]))
                .chromosomeBand(values[12]);

        if (values.length == 19) {
            builder.minRegions(Integer.valueOf(values[13]))
                    .minRegionStart(Long.valueOf(values[14]))
                    .minRegionEnd(Long.valueOf(values[15]))
                    .minRegionStartSupport(SegmentSupport.valueOf(values[16]))
                    .minRegionEndSupport(SegmentSupport.valueOf(values[17]))
                    .minRegionMethod(CopyNumberMethod.valueOf(values[18]));
        } else {
            builder.minRegions(0)
                    .minRegionStart(Long.valueOf(values[1]))
                    .minRegionEnd(Long.valueOf(values[2]))
                    .minRegionMethod(CopyNumberMethod.UNKNOWN)
                    .minRegionStartSupport(SegmentSupport.NONE)
                    .minRegionEndSupport(SegmentSupport.NONE);
        }

        return builder.build();
    }
}

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
    private static final String HEADER_PREFIX = "Chr";
    private static final String EXTENSION = ".purple.gene.cnv";

    private GeneCopyNumberFile() {
    }

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
                .add("Unused")
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
                .add("NonsenseBiallelicCount")
                .add("NonsenseNonBiallelicCount")
                .add("NonsenseNonBiallelicPloidy")
                .add("SpliceBiallelicCount")
                .add("SpliceNonBiallelicCount")
                .add("SpliceNonBiallelicPloidy")
                .add("MissenseBiallelicCount")
                .add("MissenseNonBiallelicCount")
                .add("MissenseNonBiallelicPloidy")
                .add("MinMinorAllelePloidy")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final GeneCopyNumber geneCopyNumber) {
        return new StringJoiner(DELIMITER).add(String.valueOf(geneCopyNumber.chromosome()))
                .add(String.valueOf(geneCopyNumber.start()))
                .add(String.valueOf(geneCopyNumber.end()))
                .add(String.valueOf(geneCopyNumber.gene()))
                .add(FORMAT.format(geneCopyNumber.minCopyNumber()))
                .add(FORMAT.format(geneCopyNumber.maxCopyNumber()))
                .add(String.valueOf(0))
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
                .add(String.valueOf(geneCopyNumber.nonsenseBiallelicCount()))
                .add(String.valueOf(geneCopyNumber.nonsenseNonBiallelicCount()))
                .add(FORMAT.format(geneCopyNumber.nonsenseNonBiallelicPloidy()))
                .add(String.valueOf(geneCopyNumber.spliceBiallelicCount()))
                .add(String.valueOf(geneCopyNumber.spliceNonBiallelicCount()))
                .add(FORMAT.format(geneCopyNumber.spliceNonBiallelicPloidy()))
                .add(String.valueOf(geneCopyNumber.missenseBiallelicCount()))
                .add(String.valueOf(geneCopyNumber.missenseNonBiallelicCount()))
                .add(FORMAT.format(geneCopyNumber.missenseNonBiallelicPloidy()))
                .add(FORMAT.format(geneCopyNumber.minMinorAllelePloidy()))
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
                .somaticRegions(Integer.valueOf(values[7]))
                .germlineHomRegions(Integer.valueOf(values[8]))
                .germlineHet2HomRegions(Integer.valueOf(values[9]))
                .transcriptID(values[10])
                .transcriptVersion(Integer.valueOf(values[11]))
                .chromosomeBand(values[12]);

        builder.minRegions(0)
                .minRegionStart(Long.valueOf(values[1]))
                .minRegionEnd(Long.valueOf(values[2]))
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEndSupport(SegmentSupport.NONE)
                .nonsenseBiallelicCount(0)
                .nonsenseNonBiallelicCount(0)
                .nonsenseNonBiallelicPloidy(0)
                .spliceBiallelicCount(0)
                .spliceNonBiallelicCount(0)
                .spliceNonBiallelicPloidy(0)
                .missenseBiallelicCount(0)
                .missenseNonBiallelicCount(0)
                .missenseNonBiallelicPloidy(0)
                .minMinorAllelePloidy(0);

        if (values.length >= 19) {
            builder.minRegions(Integer.valueOf(values[13]))
                    .minRegionStart(Long.valueOf(values[14]))
                    .minRegionEnd(Long.valueOf(values[15]))
                    .minRegionStartSupport(SegmentSupport.valueOf(values[16]))
                    .minRegionEndSupport(SegmentSupport.valueOf(values[17]))
                    .minRegionMethod(CopyNumberMethod.valueOf(values[18]));
        }
        if (values.length == 29) {
            builder.nonsenseBiallelicCount(Integer.valueOf(values[19]))
                    .nonsenseNonBiallelicCount(Integer.valueOf(values[20]))
                    .nonsenseNonBiallelicPloidy(Double.valueOf(values[21]))
                    .spliceBiallelicCount(Integer.valueOf(values[22]))
                    .spliceNonBiallelicCount(Integer.valueOf(values[23]))
                    .spliceNonBiallelicPloidy(Double.valueOf(values[24]))
                    .missenseBiallelicCount(Integer.valueOf(values[25]))
                    .missenseNonBiallelicCount(Integer.valueOf(values[26]))
                    .missenseNonBiallelicPloidy(Double.valueOf(values[27]))
                    .minMinorAllelePloidy(Double.valueOf(values[28]));
        }

        return builder.build();
    }
}

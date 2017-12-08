package com.hartwig.hmftools.common.gene;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

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
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final GeneCopyNumber ratio) {
        return new StringJoiner(DELIMITER).add(String.valueOf(ratio.chromosome()))
                .add(String.valueOf(ratio.start()))
                .add(String.valueOf(ratio.end()))
                .add(String.valueOf(ratio.gene()))
                .add(String.valueOf(ratio.minCopyNumber()))
                .add(String.valueOf(ratio.maxCopyNumber()))
                .add(String.valueOf(ratio.meanCopyNumber()))
                .add(String.valueOf(ratio.somaticRegions()))
                .add(String.valueOf(ratio.germlineHomRegions()))
                .add(String.valueOf(ratio.germlineHet2HomRegions()))
                .add(String.valueOf(ratio.transcriptID()))
                .add(String.valueOf(ratio.transcriptVersion()))
                .add(String.valueOf(ratio.chromosomeBand()))
                .toString();
    }

    @NotNull
    @VisibleForTesting
    static List<GeneCopyNumber> fromLines(@NotNull List<String> lines) {
        return lines.stream().filter(x -> !x.startsWith(HEADER_PREFIX)).map(GeneCopyNumberFile::fromString).collect(toList());
    }

    @NotNull
    private static GeneCopyNumber fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableGeneCopyNumber.builder()
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
                .chromosomeBand(values[12])
                .build();
    }
}

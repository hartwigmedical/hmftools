package com.hartwig.hmftools.common.fusion;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class ReportableGeneFusionFile {

    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("0.0000");
    private static final String DELIMITER = "\t";

    private static final String FILE_EXTENSION = ".linx.reported_fusion.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<ReportableGeneFusion> read(final String filePath) throws IOException {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<ReportableGeneFusion> fusions) throws IOException {
        Files.write(new File(filename).toPath(), toLines(fusions));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<ReportableGeneFusion> fusions) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        fusions.stream().map(ReportableGeneFusionFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static List<ReportableGeneFusion> fromLines(@NotNull List<String> lines) {
        return lines.stream().filter(x -> !x.startsWith("geneStart")).map(ReportableGeneFusionFile::fromString).collect(toList());
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER).add("geneStart")
                .add("geneContextStart")
                .add("geneTranscriptStart")
                .add("geneEnd")
                .add("geneContextEnd")
                .add("geneTranscriptEnd")
                .add("ploidy")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull ReportableGeneFusion fusion) {
        return new StringJoiner(DELIMITER).add(fusion.geneStart())
                .add(fusion.geneContextStart())
                .add(fusion.geneTranscriptStart())
                .add(fusion.geneEnd())
                .add(fusion.geneContextEnd())
                .add(fusion.geneTranscriptEnd())
                .add(DECIMAL_FORMAT.format(fusion.ploidy()))
                .toString();
    }

    @NotNull
    private static ReportableGeneFusion fromString(@NotNull String line) {
        String[] values = line.split(DELIMITER);

        int index = 0;

        return ImmutableReportableGeneFusion.builder()
                .geneStart(values[index++])
                .geneContextStart(values[index++])
                .geneTranscriptStart(values[index++])
                .geneEnd(values[index++])
                .geneContextEnd(values[index++])
                .geneTranscriptEnd(values[index++])
                .ploidy(Double.valueOf(values[index++]))
                .build();
    }

    @NotNull
    public static String context(@NotNull Transcript transcript, int fusedExon) {
        switch (transcript.regionType()) {
            case UPSTREAM:
                return "Promoter Region";
            case IG:
                return "IG";
            case EXONIC:
            case INTRONIC:
                return String.format("Exon %d", fusedExon);
        }

        return String.format("ERROR: %s", transcript.regionType());
    }

    public static double fusionPloidy(double downstreamPloidy, double upstreamPloidy) {
        return (upstreamPloidy + downstreamPloidy) * 0.5;
    }
}

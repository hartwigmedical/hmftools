package com.hartwig.hmftools.common.variant.structural.annotation;

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

    public static final String FILE_EXTENSION = ".linx.fusions.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<ReportableGeneFusion> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<ReportableGeneFusion> fusions) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(fusions));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<ReportableGeneFusion> fusions)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        fusions.stream().map(ReportableGeneFusionFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static List<ReportableGeneFusion> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("geneStart")).map(ReportableGeneFusionFile::fromString).collect(toList());
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("geneStart")
                .add("geneContextStart")
                .add("geneTranscriptStart")
                .add("geneEnd")
                .add("geneContextEnd")
                .add("geneTranscriptEnd")
                .add("ploidy")
                .add("source")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull ReportableGeneFusion fusion) {
        return new StringJoiner(DELIMITER).add(String.valueOf(fusion.geneStart()))
                .add(String.valueOf(fusion.geneContextStart()))
                .add(String.valueOf(fusion.geneTranscriptStart()))
                .add(String.valueOf(fusion.geneEnd()))
                .add(String.valueOf(fusion.geneContextEnd()))
                .add(String.valueOf(fusion.geneTranscriptEnd()))
                .add(DECIMAL_FORMAT.format(fusion.ploidy()))
                .add(String.valueOf(fusion.source()))
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
                .source(values[index++])
                .build();
    }

    public static String context(final String regionType, int exon, boolean isEnd)
    {
        switch (regionType)
        {
            case "Upstream":
                return "Promoter Region";
            case "Exonic":
                return String.format("Exon %d", exon);
            case "Intronic":
                return String.format("Intron %d", isEnd ? exon - 1 : exon);
        }

        return String.format("ERROR: %s", regionType);
    }

    public static double fusionPloidy(double downstreamPloidy, double upstreamPloidy)
    {
        return (upstreamPloidy + downstreamPloidy) * 0.5;
    }
}

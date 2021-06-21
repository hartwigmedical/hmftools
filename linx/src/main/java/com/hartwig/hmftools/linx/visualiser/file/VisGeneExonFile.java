package com.hartwig.hmftools.linx.visualiser.file;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.linx.visualiser.file.VisCopyNumberFile.DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class VisGeneExonFile
{
    public final String SampleId;
    public final int ClusterId;
    public final String Gene;
    public final String Transcript;
    public final String Chromosome;
    public final String AnnotationType;
    public final int ExonRank;
    public final int ExonStart;
    public final int ExonEnd;

    public VisGeneExonFile(final String sampleId, int clusterId, final String gene, final String transcript, final String chromosome,
            final String type, int exonRank, int exonStart, int exonEnd)
    {
        SampleId = sampleId;
        ClusterId = clusterId;
        Gene = gene;
        Transcript = transcript;
        Chromosome = chromosome;
        AnnotationType = type;
        ExonRank = exonRank;
        ExonStart = exonStart;
        ExonEnd = exonEnd;
    }

    private static final String FILE_EXTENSION = ".linx.vis_gene_exon.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<VisGeneExonFile> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<VisGeneExonFile> cnDataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(cnDataList));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<VisGeneExonFile> cnDataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        cnDataList.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    static List<VisGeneExonFile> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("SampleId")).map(VisGeneExonFile::fromString).collect(toList());
    }

    @NotNull
    public static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("SampleId")
                .add("ClusterId")
                .add("Gene")
                .add("Transcript")
                .add("Chromosome")
                .add("AnnotationType")
                .add("ExonRank")
                .add("ExonStart")
                .add("ExonEnd")
                .toString();
    }

    @NotNull
    public static String toString(@NotNull final VisGeneExonFile geData)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(geData.SampleId))
                .add(String.valueOf(geData.ClusterId))
                .add(String.valueOf(geData.Gene))
                .add(String.valueOf(geData.Transcript))
                .add(String.valueOf(geData.Chromosome))
                .add(String.valueOf(geData.AnnotationType))
                .add(String.valueOf(geData.ExonRank))
                .add(String.valueOf(geData.ExonStart))
                .add(String.valueOf(geData.ExonEnd))
                .toString();
    }

    @NotNull
    private static VisGeneExonFile fromString(@NotNull final String tiData)
    {
        String[] values = tiData.split(DELIMITER);

        int index = 0;

        return new VisGeneExonFile(
                values[index++],
                Integer.parseInt(values[index++]),
                values[index++],
                values[index++],
                values[index++],
                values[index++],
                Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]));
    }

}

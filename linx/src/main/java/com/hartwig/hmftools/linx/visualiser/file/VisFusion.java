package com.hartwig.hmftools.linx.visualiser.file;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.linx.visualiser.file.VisCopyNumber.DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class VisFusion
{
    public final String SampleId;
    public final boolean Reportable;
    public final int ClusterId;

    public final String GeneNameUp;
    public final String TranscriptUp;
    public final String ChrUp;
    public final int PosUp;
    public final int StrandUp;
    public final String RegionTypeUp;
    public final int FusedExonUp;

    public final String GeneNameDown;
    public final String TranscriptDown;
    public final String ChrDown;
    public final int PosDown;
    public final int StrandDown;
    public final String RegionTypeDown;
    public final int FusedExonDown;

    public VisFusion(final String sampleId, int clusterId, boolean reportable,
            final String geneNameUp, final String transcriptUp,
            final String chrUp, int posUp, int strandUp, final String regionTypeUp, int fusedExonUp,
            final String geneNameDown, final String transcriptDown,
            final String chrDown, int posDown, int strandDown, final String regionTypeDown, int fusedExonDown)
    {
        SampleId = sampleId;
        ClusterId = clusterId;
        Reportable = reportable;

        GeneNameUp = geneNameUp;
        TranscriptUp = transcriptUp;
        ChrUp = chrUp;
        PosUp = posUp;
        StrandUp = strandUp;
        RegionTypeUp = regionTypeUp;
        FusedExonUp = fusedExonUp;

        GeneNameDown = geneNameDown;
        TranscriptDown = transcriptDown;
        ChrDown = chrDown;
        PosDown = posDown;
        StrandDown = strandDown;
        RegionTypeDown = regionTypeDown;
        FusedExonDown = fusedExonDown;
    }

    public String name()
    {
        return GeneNameUp + "_" + GeneNameDown;
    }

    private static final String FILE_EXTENSION = ".linx.vis_fusion.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<VisFusion> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<VisFusion> dataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(dataList));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<VisFusion> dataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        dataList.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    static List<VisFusion> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("SampleId")).map(VisFusion::fromString).collect(toList());
    }

    @NotNull
    public static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("SampleId")
                .add("ClusterId")
                .add("Reportable")
                .add("GeneNameUp")
                .add("TranscriptUp")
                .add("ChrUp")
                .add("PosUp")
                .add("StrandUp")
                .add("RegionTypeUp")
                .add("FusedExonUp")
                .add("GeneNameDown")
                .add("TranscriptDown")
                .add("ChrDown")
                .add("PosDown")
                .add("StrandDown")
                .add("RegionTypeDown")
                .add("FusedExonDown")
                .toString();
    }

    @NotNull
    public static String toString(@NotNull final VisFusion data)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(data.SampleId))
                .add(String.valueOf(data.ClusterId))
                .add(String.valueOf(data.Reportable))
                .add(String.valueOf(data.GeneNameUp))
                .add(String.valueOf(data.TranscriptUp))
                .add(String.valueOf(data.ChrUp))
                .add(String.valueOf(data.PosUp))
                .add(String.valueOf(data.StrandUp))
                .add(String.valueOf(data.RegionTypeUp))
                .add(String.valueOf(data.FusedExonUp))
                .add(String.valueOf(data.GeneNameDown))
                .add(String.valueOf(data.TranscriptDown))
                .add(String.valueOf(data.ChrDown))
                .add(String.valueOf(data.PosDown))
                .add(String.valueOf(data.StrandDown))
                .add(String.valueOf(data.RegionTypeDown))
                .add(String.valueOf(data.FusedExonDown))
                .toString();
    }

    @NotNull
    private static VisFusion fromString(@NotNull final String data)
    {
        String[] values = data.split(DELIMITER);

        int index = 0;

        return new VisFusion(
                values[index++],
                Integer.parseInt(values[index++]),
                Boolean.parseBoolean(values[index++]),
                values[index++],
                values[index++],
                values[index++],
                Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]),
                values[index++],
                Integer.parseInt(values[index++]),
                values[index++],
                values[index++],
                values[index++],
                Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]),
                values[index++],
                Integer.parseInt(values[index++]));
    }

}

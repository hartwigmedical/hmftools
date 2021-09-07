package com.hartwig.hmftools.common.sv.linx;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class LinxViralInsertion
{
    public final String SampleId;
    public final int SvId;
    public final String VirusId;
    public final String VirusName;

    public LinxViralInsertion(final String sampleId, int svId, final String virusId, final String virusName)
    {
        SampleId = sampleId;
        SvId = svId;
        VirusId = virusId;
        VirusName = virusName;
    }

    private static final String FILE_EXTENSION = ".linx.viral_inserts.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<LinxViralInsertion> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<LinxViralInsertion> inserts) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(inserts));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<LinxViralInsertion> viruses)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        viruses.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static List<LinxViralInsertion> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("SampleId")).map(LinxViralInsertion::fromString).collect(toList());
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(LinxCluster.DELIMITER)
                .add("SampleId")
                .add("SvId")
                .add("VirusId")
                .add("VirusName")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final LinxViralInsertion virus)
    {
        return new StringJoiner(LinxCluster.DELIMITER)
                .add(String.valueOf(virus.SampleId))
                .add(String.valueOf(virus.SvId))
                .add(String.valueOf(virus.VirusId))
                .add(String.valueOf(virus.VirusName))
                .toString();
    }

    @NotNull
    private static LinxViralInsertion fromString(@NotNull final String virus)
    {
        String[] values = virus.split(LinxCluster.DELIMITER);

        int index = 0;

        return new LinxViralInsertion(
                values[index++],
                Integer.parseInt(values[index++]),
                values[index++],
                values[index++]);
    }
}

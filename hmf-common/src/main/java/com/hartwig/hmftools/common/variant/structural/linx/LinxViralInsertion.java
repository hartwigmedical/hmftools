package com.hartwig.hmftools.common.variant.structural.linx;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.variant.structural.linx.LinxClusterFile.DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class LinxViralInsertFile
{
    public final String SampleId;
    public final int SvId;
    public final String VirusId;
    public final String VirusName;

    public LinxViralInsertFile(final String sampleId, int svId, final String virusId, final String virusName)
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
    public static List<LinxViralInsertFile> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<LinxViralInsertFile> inserts) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(inserts));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<LinxViralInsertFile> viruses)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        viruses.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static List<LinxViralInsertFile> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("SampleId")).map(LinxViralInsertFile::fromString).collect(toList());
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("SampleId")
                .add("SvId")
                .add("VirusId")
                .add("VirusName")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final LinxViralInsertFile virus)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(virus.SampleId))
                .add(String.valueOf(virus.SvId))
                .add(String.valueOf(virus.VirusId))
                .add(String.valueOf(virus.VirusName))
                .toString();
    }

    @NotNull
    private static LinxViralInsertFile fromString(@NotNull final String virus)
    {
        String[] values = virus.split(DELIMITER);

        int index = 0;

        return new LinxViralInsertFile(
                values[index++],
                Integer.parseInt(values[index++]),
                values[index++],
                values[index++]);
    }
}

package com.hartwig.hmftools.common.linx;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

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

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static List<LinxViralInsertion> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(final String filename, List<LinxViralInsertion> inserts) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(inserts));
    }

    private static List<String> toLines(final List<LinxViralInsertion> viruses)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        viruses.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    private static List<LinxViralInsertion> fromLines(List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("SampleId")).map(LinxViralInsertion::fromString).collect(toList());
    }

    private static String header()
    {
        return new StringJoiner(LinxCluster.DELIMITER)
                .add("SampleId")
                .add("SvId")
                .add("VirusId")
                .add("VirusName")
                .toString();
    }

    private static String toString(final LinxViralInsertion virus)
    {
        return new StringJoiner(LinxCluster.DELIMITER)
                .add(String.valueOf(virus.SampleId))
                .add(String.valueOf(virus.SvId))
                .add(String.valueOf(virus.VirusId))
                .add(String.valueOf(virus.VirusName))
                .toString();
    }

    private static LinxViralInsertion fromString(final String virus)
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

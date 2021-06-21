package com.hartwig.hmftools.linx.visualiser.file;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class VisCopyNumberFile
{
    public final String SampleId;
    public final String Chromosome;
    public final int Start;
    public final int End;
    public final double CopyNumber;
    public final double BAF;

    public VisCopyNumberFile(final String sampleId, final String chromosome, int start, int end, double copyNumber, double baf)
    {
        SampleId = sampleId;
        Chromosome = chromosome;
        Start = start;
        End = end;
        CopyNumber = copyNumber;
        BAF = baf;
    }

    public static final String DELIMITER = "\t";
    public static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final String FILE_EXTENSION = ".linx.vis_copy_number.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<VisCopyNumberFile> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<VisCopyNumberFile> cnDataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(cnDataList));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<VisCopyNumberFile> cnDataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        cnDataList.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    static List<VisCopyNumberFile> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("SampleId")).map(VisCopyNumberFile::fromString).collect(toList());
    }

    @NotNull
    public static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("SampleId")
                .add("Chromosome")
                .add("Start")
                .add("End")
                .add("CopyNumber")
                .add("BAF")
                .toString();
    }

    @NotNull
    public static String toString(@NotNull final VisCopyNumberFile cnData)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(cnData.SampleId))
                .add(String.valueOf(cnData.Chromosome))
                .add(String.valueOf(cnData.Start))
                .add(String.valueOf(cnData.End))
                .add(FORMAT.format(cnData.CopyNumber))
                .add(FORMAT.format(cnData.BAF))
                .toString();
    }

    @NotNull
    private static VisCopyNumberFile fromString(@NotNull final String tiData)
    {
        String[] values = tiData.split(DELIMITER);

        int index = 0;

        return new VisCopyNumberFile(
                values[index++],
                values[index++],
                Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]),
                Double.parseDouble(values[index++]),
                Double.parseDouble(values[index++]));
    }
}

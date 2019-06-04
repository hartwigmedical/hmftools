package com.hartwig.hmftools.svanalysis.visual;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.svanalysis.visual.VisCopyNumberFile.DELIMITER;
import static com.hartwig.hmftools.svanalysis.visual.VisCopyNumberFile.HEADER_PREFIX;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class VisSegmentFile
{
    public final String SampleId;
    public final int ClusterId;
    public final int ChainId;
    public final String Chromosome;
    public final String PosStart;
    public final String PosEnd;
    public final int TraverseCount;


    public VisSegmentFile(final String sampleId, int clusterId, int chainId, final String chromosome,
            final String posStart, final String posEnd, int traverseCount)
    {
        SampleId = sampleId;
        ClusterId = clusterId;
        ChainId = chainId;
        Chromosome = chromosome;
        PosStart = posStart;
        PosEnd = posEnd;
        TraverseCount = traverseCount;
    }

    private static final String FILE_EXTENSION = ".linx.vis_segments.csv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<VisSegmentFile> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<VisSegmentFile> svDataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(svDataList));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<VisSegmentFile> segments)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        segments.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    static List<VisSegmentFile> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith(HEADER_PREFIX)).map(VisSegmentFile::fromString).collect(toList());
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER, HEADER_PREFIX,"")
                .add("SampleId")
                .add("ClusterId")
                .add("ChainId")
                .add("Chromosome")
                .add("PosStart")
                .add("PosEnd")
                .add("TraverseCount")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final VisSegmentFile segment)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(segment.SampleId))
                .add(String.valueOf(segment.ClusterId))
                .add(String.valueOf(segment.ChainId))
                .add(String.valueOf(segment.Chromosome))
                .add(String.valueOf(segment.PosStart))
                .add(String.valueOf(segment.PosEnd))
                .add(String.valueOf(segment.TraverseCount))
                .toString();
    }

    @NotNull
    private static VisSegmentFile fromString(@NotNull final String segment)
    {
        String[] values = segment.split(DELIMITER);

        int index = 0;

        return new VisSegmentFile(
                values[index++],
                Integer.valueOf(values[index++]),
                Integer.valueOf(values[index++]),
                values[index++],
                values[index++],
                values[index++],
                Integer.valueOf(values[index++]));
    }

}

package com.hartwig.hmftools.linx.visualiser.file;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.linx.visualiser.file.VisCopyNumber.DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.linx.visualiser.circos.SegmentTerminal;

import org.jetbrains.annotations.NotNull;

public class VisSegment implements GenomeRegion
{
    public final String SampleId;
    public final int ClusterId;
    public final int ChainId;
    public final String Chromosome;
    public final String PosStart;
    public final String PosEnd;
    public final double LinkPloidy;
    public final boolean InDoubleMinute;

    public int Frame;
    public int Track;

    private int mStart;
    private int mEnd;
    private SegmentTerminal mTerminalStart;
    private SegmentTerminal mTerminalEnd;

    public VisSegment(
            final String sampleId, int clusterId, int chainId, final String chromosome, final String posStart, final String posEnd,
            double linkPloidy, boolean inDM)
    {
        SampleId = sampleId;
        ClusterId = clusterId;
        ChainId = chainId;
        Chromosome = chromosome;
        PosStart = posStart;
        PosEnd = posEnd;
        LinkPloidy = linkPloidy;
        InDoubleMinute = inDM;
        Frame = 0;
        Track = 0;

        mTerminalStart = SegmentTerminal.fromString(PosStart);
        mTerminalEnd = SegmentTerminal.fromString(PosEnd);
        mStart = SegmentTerminal.fromString(PosStart) == SegmentTerminal.NONE ? Integer.valueOf(PosStart) : Integer.valueOf(PosEnd);
        mEnd = SegmentTerminal.fromString(PosEnd) == SegmentTerminal.NONE ? Integer.valueOf(PosEnd) : Integer.valueOf(PosStart);
    }

    @Override
    public String chromosome() { return Chromosome; }

    @Override
    public int start() { return mStart; }

    @Override
    public int end() { return mEnd; }

    public void setStart(int pos) { mStart = pos; }
    public void setEnd(int pos) { mEnd = pos; }

    public SegmentTerminal startTerminal() { return mTerminalStart; }
    public SegmentTerminal endTerminal() { return mTerminalEnd; };

    public void setTerminalStart(final SegmentTerminal terminal) { mTerminalStart = terminal; }
    public void setTerminalEnd(final SegmentTerminal terminal) { mTerminalEnd = terminal; }

    public static VisSegment from(final VisSegment other)
    {
        VisSegment newSegment = new VisSegment(
                other.SampleId, other.ClusterId, other.ChainId, other.Chromosome, other.PosStart, other.PosEnd,
                other.LinkPloidy, other.InDoubleMinute);

        newSegment.setTerminalStart(other.startTerminal());
        newSegment.setTerminalEnd(other.endTerminal());
        newSegment.Track = other.Track;
        newSegment.Frame = other.Frame;
        newSegment.setStart(other.start());
        newSegment.setEnd(other.end());
        newSegment.setTerminalStart(other.startTerminal());
        newSegment.setTerminalEnd(other.endTerminal());
        return newSegment;
    }

    private static final String FILE_EXTENSION = ".linx.vis_segments.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<VisSegment> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<VisSegment> svDataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(svDataList));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<VisSegment> segments)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        segments.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    static List<VisSegment> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("SampleId")).map(VisSegment::fromString).collect(toList());
    }

    @NotNull
    public static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("SampleId")
                .add("ClusterId")
                .add("ChainId")
                .add("Chromosome")
                .add("PosStart")
                .add("PosEnd")
                .add("LinkPloidy")
                .add("InDoubleMinute")
                .toString();
    }

    @NotNull
    public static String toString(@NotNull final VisSegment segment)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(segment.SampleId))
                .add(String.valueOf(segment.ClusterId))
                .add(String.valueOf(segment.ChainId))
                .add(String.valueOf(segment.Chromosome))
                .add(String.valueOf(segment.PosStart))
                .add(String.valueOf(segment.PosEnd))
                .add(String.format("%.4f",segment.LinkPloidy))
                .add(String.valueOf(segment.InDoubleMinute))
                .toString();
    }

    @NotNull
    private static VisSegment fromString(@NotNull final String segment)
    {
        String[] values = segment.split(DELIMITER);

        int index = 0;

        return new VisSegment(
                values[index++],
                Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]),
                values[index++],
                values[index++],
                values[index++],
                Double.parseDouble(values[index++]),
                index < values.length ? Boolean.parseBoolean(values[index++]) : false);
    }

}

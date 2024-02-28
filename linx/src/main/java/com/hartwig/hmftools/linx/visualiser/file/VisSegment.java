package com.hartwig.hmftools.linx.visualiser.file;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getBoolValue;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getDoubleValue;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getIntValue;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getValue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.linx.visualiser.circos.SegmentTerminal;

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
    private static final String GERMLINE_FILE_EXTENSION = ".linx.germline.vis_segments.tsv";

    public static String generateFilename(final String basePath, final String sample, boolean isGermline)
    {
        return basePath + File.separator + sample + (isGermline ? GERMLINE_FILE_EXTENSION : FILE_EXTENSION);
    }

    public static List<VisSegment> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(final String filename, List<VisSegment> svDataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(svDataList));
    }

    private static List<String> toLines(final List<VisSegment> segments)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        segments.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    private static List<VisSegment> fromLines(final List<String> lines)
    {
        String header = lines.get(0);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        List<VisSegment> data = Lists.newArrayList();

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);

            data.add(new VisSegment(
                    getValue(fieldsIndexMap, FLD_SAMPLE_ID, "", values),
                    getIntValue(fieldsIndexMap, "ClusterId", 0, values),
                    getIntValue(fieldsIndexMap, "ChainId", 0, values),
                    getValue(fieldsIndexMap, "Chromosome", "", values),
                    getValue(fieldsIndexMap, "PosStart", "", values),
                    getValue(fieldsIndexMap, "PosEnd", "", values),
                    getDoubleValue(fieldsIndexMap, "LinkPloidy", 0, values),
                    getBoolValue(fieldsIndexMap, "InDoubleMinute", false, values)));

        }

        return data;
    }

    public static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add("ClusterId")
                .add("ChainId")
                .add("Chromosome")
                .add("PosStart")
                .add("PosEnd")
                .add("LinkPloidy")
                .add("InDoubleMinute")
                .toString();
    }

    public static String toString(final VisSegment segment)
    {
        return new StringJoiner(TSV_DELIM)
                .add(String.valueOf(segment.ClusterId))
                .add(String.valueOf(segment.ChainId))
                .add(String.valueOf(segment.Chromosome))
                .add(String.valueOf(segment.PosStart))
                .add(String.valueOf(segment.PosEnd))
                .add(String.format("%.4f",segment.LinkPloidy))
                .add(String.valueOf(segment.InDoubleMinute))
                .toString();
    }
}

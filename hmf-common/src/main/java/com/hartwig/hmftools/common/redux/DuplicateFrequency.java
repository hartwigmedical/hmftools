package com.hartwig.hmftools.common.redux;

import static com.hartwig.hmftools.common.metrics.BamMetricsCommon.COUNT_COLUMN;
import static com.hartwig.hmftools.common.redux.ReduxCommon.REDUX_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public class DuplicateFrequency
{
    public final int ReadCount;
    public long Count;
    public long DualStrandCount;

    public DuplicateFrequency(final int readCount)
    {
        this(readCount, 0, NO_DUAL_STRAND);
    }

    public DuplicateFrequency(final int readCount, final long count)
    {
        this(readCount, count, NO_DUAL_STRAND);
    }

    public DuplicateFrequency(final int readCount, final long count, final long dualStrandCount)
    {
        ReadCount = readCount;
        Count = count;
        DualStrandCount = dualStrandCount;
    }

    public boolean hasDualStrandCount() { return DualStrandCount >= NO_DUAL_STRAND; }

    private static final String FILE_EXTENSION = ".duplicate_freq.tsv";
    public static final int NO_DUAL_STRAND = -1;

    public static String generateFilename(final String basePath, final String sampleId)
    {
        return checkAddDirSeparator(basePath) + sampleId + REDUX_FILE_ID + FILE_EXTENSION;
    }

    public static final String DUP_READ_COUNT_COLUMN = "DuplicateReadCount";
    public static final String DUAL_STRAND_COUNT_COLUMN = "DualStrandFrequency";

    public static void write(final String filename, final List<DuplicateFrequency> frequencies, boolean includeDualStrand) throws IOException
    {
        BufferedWriter writer = createBufferedWriter(filename, false);

        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(DUP_READ_COUNT_COLUMN);
        sj.add(COUNT_COLUMN);

        if(includeDualStrand)
            sj.add(DUAL_STRAND_COUNT_COLUMN);

        writer.write(sj.toString());
        writer.newLine();

        for(DuplicateFrequency frequency : frequencies)
        {
            sj = new StringJoiner(TSV_DELIM);
            sj.add(String.valueOf(frequency.ReadCount));
            sj.add(String.valueOf(frequency.Count));

            if(includeDualStrand)
                sj.add(String.valueOf(frequency.DualStrandCount));

            writer.write(sj.toString());
            writer.newLine();
        }

        writer.close();
    }

    public static List<DuplicateFrequency> read(final String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());

        String header = lines.get(0);
        lines.remove(0);

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

        // parse dynamic coverage depth percentages
        List<DuplicateFrequency> frequencies = Lists.newArrayListWithCapacity(lines.size());

        int dupReadIndex = fieldsIndexMap.get(DUP_READ_COUNT_COLUMN);
        int countIndex = fieldsIndexMap.get(COUNT_COLUMN);
        Integer dualStrandIndex = fieldsIndexMap.get(DUAL_STRAND_COUNT_COLUMN);

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            int dupReads = Integer.parseInt(values[dupReadIndex]);
            long count = Long.parseLong(values[countIndex]);

            long dualStrand = dualStrandIndex != null ? Long.parseLong(values[dualStrandIndex]) : NO_DUAL_STRAND;

            frequencies.add(new DuplicateFrequency(dupReads, count, dualStrand));
        }

        return frequencies;
    }
}

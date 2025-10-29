package com.hartwig.hmftools.common.metrics;

import static com.hartwig.hmftools.common.metrics.BamMetricsCommon.BAM_METRICS_FILE_ID;
import static com.hartwig.hmftools.common.metrics.BamMetricsCommon.COUNT_COLUMN;
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

public class BamMetricFragmentLength
{
    public final List<ValueFrequency> FragmentLengths;

    private static final String FILE_EXTENSION = ".frag_length.tsv";

    public static String generateFilename(final String basePath, final String sampleId)
    {
        return checkAddDirSeparator(basePath) + sampleId + BAM_METRICS_FILE_ID + FILE_EXTENSION;
    }

    public BamMetricFragmentLength(final List<ValueFrequency> fragmentLengths)
    {
        FragmentLengths = fragmentLengths;
    }

    public static final String FRAG_LENGTH_COLUMN = "FragmentLength";

    public static void write(final String filename, final List<ValueFrequency> fragmentLengths) throws IOException
    {
        BufferedWriter writer = createBufferedWriter(filename, false);

        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(FRAG_LENGTH_COLUMN).add(COUNT_COLUMN);
        writer.write(sj.toString());
        writer.newLine();

        for(ValueFrequency lengthCount : fragmentLengths)
        {
            sj = new StringJoiner(TSV_DELIM);
            sj.add(String.valueOf(lengthCount.Value)).add(String.valueOf(lengthCount.Count));
            writer.write(sj.toString());
            writer.newLine();
        }

        writer.close();
    }

    public static BamMetricFragmentLength read(final String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());

        String header = lines.get(0);
        lines.remove(0);

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

        // parse dynamic coverage depth percentages
        List<ValueFrequency> fragmentLengths = Lists.newArrayListWithCapacity(lines.size());

        int fragLengthIndex = fieldsIndexMap.get(FRAG_LENGTH_COLUMN);
        int countIndex = fieldsIndexMap.get(COUNT_COLUMN);

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            int fragLength = Integer.parseInt(values[fragLengthIndex]);
            long count = Long.parseLong(values[countIndex]);

            fragmentLengths.add(new ValueFrequency(fragLength, count));
        }

        return new BamMetricFragmentLength(fragmentLengths);
    }
}

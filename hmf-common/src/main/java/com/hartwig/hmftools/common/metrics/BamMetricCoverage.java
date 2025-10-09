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

public class BamMetricCoverage
{
    public final List<ValueFrequency> Coverage;

    private static final String FILE_EXTENSION = ".coverage.tsv";

    public static String generateFilename(final String basePath, final String sampleId)
    {
        return checkAddDirSeparator(basePath) + sampleId + BAM_METRICS_FILE_ID + FILE_EXTENSION;
    }

    public BamMetricCoverage(final List<ValueFrequency> coverageLevels)
    {
        Coverage = coverageLevels;
    }

    public static final String COVERAGE_COLUMN = "Coverage";

    public static void write(final String filename, final List<ValueFrequency> coverageLevels) throws IOException
    {
        BufferedWriter writer = createBufferedWriter(filename, false);

        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(COVERAGE_COLUMN).add(COUNT_COLUMN);
        writer.write(sj.toString());
        writer.newLine();

        for(ValueFrequency coverageFrequency : coverageLevels)
        {
            sj = new StringJoiner(TSV_DELIM);
            sj.add(String.valueOf(coverageFrequency.Value)).add(String.valueOf(coverageFrequency.Count));
            writer.write(sj.toString());
            writer.newLine();
        }

        writer.close();
    }

    public static BamMetricCoverage read(final String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());

        String header = lines.get(0);
        lines.remove(0);

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

        // parse dynamic coverage depth percentages
        List<ValueFrequency> coverageLevels = Lists.newArrayListWithCapacity(lines.size());

        int coverageIndex = fieldsIndexMap.get(COVERAGE_COLUMN);
        int countIndex = fieldsIndexMap.get(COUNT_COLUMN);

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            int coverage = Integer.parseInt(values[coverageIndex]);
            long count = Long.parseLong(values[countIndex]);

            coverageLevels.add(new ValueFrequency(coverage, count));
        }

        return new BamMetricCoverage(coverageLevels);
    }
}

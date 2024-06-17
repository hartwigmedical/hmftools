package com.hartwig.hmftools.common.metrics;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import org.jetbrains.annotations.NotNull;

// Picard WgsMetrics file output, internally superceded with BamMetricsSummary

public final class WGSMetricsFile
{
    public static final String GENOME_TERRITORY_COLUMN = "GENOME_TERRITORY";
    public static final String MEAN_COVERAGE_COLUMN = "MEAN_COVERAGE";
    public static final String SD_COVERAGE_COLUMN = "SD_COVERAGE";
    public static final String MEDIAN_COVERAGE_COLUMN = "MEDIAN_COVERAGE";
    public static final String MAD_COVERAGE_COLUMN = "MAD_COVERAGE";
    public static final String PCT_EXC_ADAPTER_COLUMN = "PCT_EXC_ADAPTER";
    public static final String PCT_EXC_MAPQ_COLUMN = "PCT_EXC_MAPQ";
    public static final String PCT_EXC_DUPE_COLUMN = "PCT_EXC_DUPE";
    public static final String PCT_EXC_UNPAIRED_COLUMN = "PCT_EXC_UNPAIRED";
    public static final String PCT_EXC_BASEQ_COLUMN = "PCT_EXC_BASEQ";
    public static final String PCT_EXC_OVERLAP_COLUMN = "PCT_EXC_OVERLAP";
    public static final String PCT_EXC_CAPPED_COLUMN = "PCT_EXC_CAPPED";
    public static final String PCT_EXC_TOTAL_COLUMN = "PCT_EXC_TOTAL";

    private static final String COVERAGE_1X_COLUMN = "PCT_1X";
    private static final String COVERAGE_10X_COLUMN = "PCT_10X";
    private static final String COVERAGE_20X_COLUMN = "PCT_20X";
    private static final String COVERAGE_30X_COLUMN = "PCT_30X";
    private static final String COVERAGE_60X_COLUMN = "PCT_60X";

    public static final String FILE_EXTENSION = ".wgsmetrics";

    public static String generateFilename(final String basePath, final String sampleId)
    {
        return checkAddDirSeparator(basePath) + sampleId + FILE_EXTENSION;
    }

    @NotNull
    public static WGSMetrics read(final String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());

        String headerLine = null;
        String valuesLine = null;

        for(int i = 0; i < lines.size() - 1; ++i)
        {
            if(lines.get(i).startsWith(GENOME_TERRITORY_COLUMN))
            {
                headerLine = lines.get(i);
                valuesLine = lines.get(i + 1);
                break;
            }
        }

        if(headerLine == null)
        {
            throw new IOException("invalid WGS metrics file: " + filename);
        }

        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(headerLine, TSV_DELIM);
        String[] values = valuesLine.split(TSV_DELIM, -1);

        // NOTE: adapter and 1x coverage not exist in older versions

        return ImmutableWGSMetrics.builder()
                .meanCoverage(Double.parseDouble(values[fieldsIndexMap.get(MEAN_COVERAGE_COLUMN)]))
                .sdCoverage(Double.parseDouble(values[fieldsIndexMap.get(SD_COVERAGE_COLUMN)]))
                .medianCoverage((int) Double.parseDouble(values[fieldsIndexMap.get(MEDIAN_COVERAGE_COLUMN)]))
                .madCoverage((int) Double.parseDouble(values[fieldsIndexMap.get(MAD_COVERAGE_COLUMN)]))
                .pctExcAdapter(fieldsIndexMap.containsKey(PCT_EXC_ADAPTER_COLUMN) ?
                        Double.parseDouble(values[fieldsIndexMap.get(PCT_EXC_ADAPTER_COLUMN)]) : null)
                .pctExcMapQ(Double.parseDouble(values[fieldsIndexMap.get(PCT_EXC_MAPQ_COLUMN)]))
                .pctExcDupe(Double.parseDouble(values[fieldsIndexMap.get(PCT_EXC_DUPE_COLUMN)]))
                .pctExcUnpaired(Double.parseDouble(values[fieldsIndexMap.get(PCT_EXC_UNPAIRED_COLUMN)]))
                .pctExcBaseQ(Double.parseDouble(values[fieldsIndexMap.get(PCT_EXC_BASEQ_COLUMN)]))
                .pctExcOverlap(Double.parseDouble(values[fieldsIndexMap.get(PCT_EXC_OVERLAP_COLUMN)]))
                .pctExcCapped(Double.parseDouble(values[fieldsIndexMap.get(PCT_EXC_CAPPED_COLUMN)]))
                .pctExcTotal(Double.parseDouble(values[fieldsIndexMap.get(PCT_EXC_TOTAL_COLUMN)]))
                .coverage1xPercentage(fieldsIndexMap.containsKey(COVERAGE_1X_COLUMN) ?
                        Double.parseDouble(values[fieldsIndexMap.get(COVERAGE_1X_COLUMN)]) : null)
                .coverage10xPercentage(Double.parseDouble(values[fieldsIndexMap.get(COVERAGE_10X_COLUMN)]))
                .coverage20xPercentage(Double.parseDouble(values[fieldsIndexMap.get(COVERAGE_20X_COLUMN)]))
                .coverage30xPercentage(Double.parseDouble(values[fieldsIndexMap.get(COVERAGE_30X_COLUMN)]))
                .coverage60xPercentage(Double.parseDouble(values[fieldsIndexMap.get(COVERAGE_60X_COLUMN)]))
                .build();
    }
}

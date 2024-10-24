package com.hartwig.hmftools.compar.metrics;

// Picard WgsMetrics file output, internally superseded with BamMetricsSummary

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.metrics.BamMetricsSummary;
import com.hartwig.hmftools.common.metrics.ImmutableBamMetricsSummary;

import org.jetbrains.annotations.NotNull;

public final class OldWGSMetricsFile
{
    public static final String GENOME_TERRITORY_COLUMN = "GENOME_TERRITORY";
    public static final String MEAN_COVERAGE_COLUMN = "MEAN_COVERAGE";
    public static final String SD_COVERAGE_COLUMN = "SD_COVERAGE";
    public static final String MEDIAN_COVERAGE_COLUMN = "MEDIAN_COVERAGE";
    public static final String MAD_COVERAGE_COLUMN = "MAD_COVERAGE";
    public static final String PCT_EXC_MAPQ_COLUMN = "PCT_EXC_MAPQ";
    public static final String PCT_EXC_DUPE_COLUMN = "PCT_EXC_DUPE";
    public static final String PCT_EXC_UNPAIRED_COLUMN = "PCT_EXC_UNPAIRED";
    public static final String PCT_EXC_BASEQ_COLUMN = "PCT_EXC_BASEQ";
    public static final String PCT_EXC_OVERLAP_COLUMN = "PCT_EXC_OVERLAP";
    public static final String PCT_EXC_CAPPED_COLUMN = "PCT_EXC_CAPPED";

    private static final String COVERAGE_10X_COLUMN = "PCT_10X";
    private static final String COVERAGE_20X_COLUMN = "PCT_20X";
    private static final String COVERAGE_30X_COLUMN = "PCT_30X";
    private static final String COVERAGE_60X_COLUMN = "PCT_60X";

    @NotNull
    public static BamMetricsSummary read(final String filename) throws IOException
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

        List<Integer> coverageLevels = List.of(10, 20, 30, 60);
        List<Double> coveragePercents = List.of(
                Double.parseDouble(values[fieldsIndexMap.get(COVERAGE_10X_COLUMN)]),
                Double.parseDouble(values[fieldsIndexMap.get(COVERAGE_20X_COLUMN)]),
                Double.parseDouble(values[fieldsIndexMap.get(COVERAGE_30X_COLUMN)]),
                Double.parseDouble(values[fieldsIndexMap.get(COVERAGE_60X_COLUMN)])
        );
        return ImmutableBamMetricsSummary.builder()
                .totalRegionBases(-1)
                .totalReads(-1)
                .duplicateReads(-1)
                .dualStrandReads(-1)
                .meanCoverage(Double.parseDouble(values[fieldsIndexMap.get(MEAN_COVERAGE_COLUMN)]))
                .sdCoverage(Double.parseDouble(values[fieldsIndexMap.get(SD_COVERAGE_COLUMN)]))
                .medianCoverage((int) Double.parseDouble(values[fieldsIndexMap.get(MEDIAN_COVERAGE_COLUMN)]))
                .madCoverage((int) Double.parseDouble(values[fieldsIndexMap.get(MAD_COVERAGE_COLUMN)]))
                .lowMapQualPercent(Double.parseDouble(values[fieldsIndexMap.get(PCT_EXC_MAPQ_COLUMN)]))
                .duplicatePercent(Double.parseDouble(values[fieldsIndexMap.get(PCT_EXC_DUPE_COLUMN)]))
                .unpairedPercent(Double.parseDouble(values[fieldsIndexMap.get(PCT_EXC_UNPAIRED_COLUMN)]))
                .lowBaseQualPercent(Double.parseDouble(values[fieldsIndexMap.get(PCT_EXC_BASEQ_COLUMN)]))
                .overlappingReadPercent(Double.parseDouble(values[fieldsIndexMap.get(PCT_EXC_OVERLAP_COLUMN)]))
                .cappedCoveragePercent(Double.parseDouble(values[fieldsIndexMap.get(PCT_EXC_CAPPED_COLUMN)]))
                .coverageLevels(coverageLevels)
                .coveragePercents(coveragePercents)
                .build();
    }
}

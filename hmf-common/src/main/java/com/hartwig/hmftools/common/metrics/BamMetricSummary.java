package com.hartwig.hmftools.common.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.metrics.BamMetricsCommon.BAM_METRICS_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BamMetricSummary
{
    public abstract long totalRegionBases();
    public abstract long totalReads();
    public abstract long duplicateReads();
    public abstract long dualStrandReads();

    public abstract double meanCoverage();
    public abstract double sdCoverage();
    public abstract int medianCoverage();
    public abstract int madCoverage();

    public abstract double lowMapQualPercent();
    public abstract double duplicatePercent();
    public abstract double unmappedPercent();
    public abstract double lowBaseQualPercent();
    public abstract double overlappingReadPercent();
    public abstract double cappedCoveragePercent();

    public double totalFilteredPercent()
    {
        return lowMapQualPercent() + duplicatePercent() + unmappedPercent() + lowBaseQualPercent()
                + overlappingReadPercent() + cappedCoveragePercent();
    }

    public abstract List<Integer> coverageLevels();
    public abstract List<Double> coveragePercents();

    public Double getCoveragePercent(int coverageLevel)
    {
        for(int i = 0; i < coverageLevels().size(); ++i)
        {
            if(coverageLevels().get(i) == coverageLevel)
                return coveragePercents().get(i);
        }

        return null;
    }

    public double coveragePercent(int coverageLevel)
    {
        Double percent = getCoveragePercent(coverageLevel);
        return percent != null ? percent : 0;
    }

    // load coverage dynamic
    // PCT_1X	PCT_5X	PCT_10X	PCT_15X	PCT_20X	PCT_25X	PCT_30X	PCT_40X	PCT_50X	PCT_60X	PCT_70X	PCT_80X	PCT_90X	PCT_100X

    private static final String FILE_EXTENSION = ".summary.tsv";

    public static String generateFilename(final String basePath, final String sampleId)
    {
        return checkAddDirSeparator(basePath) + sampleId + BAM_METRICS_FILE_ID + FILE_EXTENSION;
    }

    public static final String TOTAL_REGION_COLUMN = "TotalRegionBases";
    public static final String TOTAL_READS_COLUMN = "TotalReads";
    public static final String DUPLICATES_COLUMN = "DuplicateReads";
    public static final String DUAL_STRAND_COLUMN = "DualStrandReads";
    public static final String MEAN_COVERAGE_COLUMN = "MeanCoverage";
    public static final String SD_COVERAGE_COLUMN = "StdDevCoverage";
    public static final String MEDIAN_COVERAGE_COLUMN = "MedianCoverage";
    public static final String MAD_COVERAGE_COLUMN = "MadCoverage";
    public static final String LOW_MAPQ_COLUMN = "LowMapQualPercent";
    public static final String DUPLICATE_PCT_COLUMN = "DuplicatePercent";
    public static final String UNMAPPED_COLUMN = "UnmappedPercent";
    public static final String LOW_BASEQ_COLUMN = "LowBaseQualPercent";
    public static final String OVERLAP_READ_COLUMN = "OverlappingReadPercent";
    public static final String CAPPED_COVERAGE_COLUMN = "CappedCoverage";

    private static final String UNPAIRED_COLUMN = "UnpairedPercent"; // decprecated in v1.5

    private static final String DEPTH_COVERAGE = "DepthCoverage";

    public void write(final String filename) throws IOException
    {
        List<String> lines = Lists.newArrayList();

        StringJoiner header = new StringJoiner(TSV_DELIM);
        header.add(TOTAL_REGION_COLUMN);
        header.add(TOTAL_READS_COLUMN);
        header.add(DUPLICATES_COLUMN);
        header.add(DUAL_STRAND_COLUMN);
        header.add(MEAN_COVERAGE_COLUMN);
        header.add(SD_COVERAGE_COLUMN);
        header.add(MEDIAN_COVERAGE_COLUMN);
        header.add(MAD_COVERAGE_COLUMN);
        header.add(LOW_MAPQ_COLUMN);
        header.add(DUPLICATE_PCT_COLUMN);
        header.add(UNMAPPED_COLUMN);
        header.add(LOW_BASEQ_COLUMN);
        header.add(OVERLAP_READ_COLUMN);
        header.add(CAPPED_COVERAGE_COLUMN);

        for(int i = 0; i < coverageLevels().size(); ++i)
        {
            header.add(format("%s_%d", DEPTH_COVERAGE, coverageLevels().get(i)));
        }

        lines.add(header.toString());

        StringJoiner values = new StringJoiner(TSV_DELIM);
        values.add(String.valueOf(totalRegionBases()));
        values.add(String.valueOf(totalReads()));
        values.add(String.valueOf(duplicateReads()));
        values.add(String.valueOf(dualStrandReads()));

        values.add(format("%.3f", meanCoverage()));
        values.add(format("%.3f", sdCoverage()));
        values.add(format("%d", medianCoverage()));
        values.add(format("%d", madCoverage()));

        values.add(format("%.5f", lowMapQualPercent()));
        values.add(format("%.5f", duplicatePercent()));
        values.add(format("%.5f", unmappedPercent()));
        values.add(format("%.5f", lowBaseQualPercent()));
        values.add(format("%.5f", overlappingReadPercent()));
        values.add(format("%.5f", cappedCoveragePercent()));

        for(int i = 0; i < coveragePercents().size(); ++i)
        {
            values.add(format("%.5f", coveragePercents().get(i)));
        }

        lines.add(values.toString());

        Files.write(new File(filename).toPath(), lines);
    }

    public static BamMetricSummary read(final String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());

        String header = lines.get(0);
        String valuesLine = lines.get(1);

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        String[] values = valuesLine.split(TSV_DELIM, -1);

        // parse dynamic coverage depth percentages
        List<Integer> coverageLevels = Lists.newArrayList();
        List<Double> coveragePercents = Lists.newArrayList();

        String[] columns = header.split(TSV_DELIM,-1);

        for(int i = 0; i < columns.length; ++i)
        {
            String column = columns[i];

            if(column.startsWith(DEPTH_COVERAGE))
            {
                int depthLevel = Integer.parseInt(column.split("_", 2)[1]);
                coverageLevels.add(depthLevel);
                coveragePercents.add(Double.parseDouble(values[i]));
            }
        }

        int unmappedIndex = fieldsIndexMap.containsKey(UNMAPPED_COLUMN) ?
                fieldsIndexMap.get(UNMAPPED_COLUMN) : fieldsIndexMap.get(UNPAIRED_COLUMN);

        return ImmutableBamMetricSummary.builder()
                .totalRegionBases(Long.parseLong(values[fieldsIndexMap.get(TOTAL_REGION_COLUMN)]))
                .totalReads(Long.parseLong(values[fieldsIndexMap.get(TOTAL_READS_COLUMN)]))
                .duplicateReads(Long.parseLong(values[fieldsIndexMap.get(DUPLICATES_COLUMN)]))
                .dualStrandReads(Long.parseLong(values[fieldsIndexMap.get(DUAL_STRAND_COLUMN)]))
                .meanCoverage(Double.parseDouble(values[fieldsIndexMap.get(MEAN_COVERAGE_COLUMN)]))
                .sdCoverage(Double.parseDouble(values[fieldsIndexMap.get(SD_COVERAGE_COLUMN)]))
                .medianCoverage((int)Double.parseDouble(values[fieldsIndexMap.get(MEDIAN_COVERAGE_COLUMN)]))
                .madCoverage((int)Double.parseDouble(values[fieldsIndexMap.get(MAD_COVERAGE_COLUMN)]))
                .lowMapQualPercent(Double.parseDouble(values[fieldsIndexMap.get(LOW_MAPQ_COLUMN)]))
                .duplicatePercent(Double.parseDouble(values[fieldsIndexMap.get(DUPLICATE_PCT_COLUMN)]))
                .unmappedPercent(Double.parseDouble(values[unmappedIndex]))
                .lowBaseQualPercent(Double.parseDouble(values[fieldsIndexMap.get(LOW_BASEQ_COLUMN)]))
                .overlappingReadPercent(Double.parseDouble(values[fieldsIndexMap.get(OVERLAP_READ_COLUMN)]))
                .cappedCoveragePercent(Double.parseDouble(values[fieldsIndexMap.get(CAPPED_COVERAGE_COLUMN)]))
                .coverageLevels(coverageLevels)
                .coveragePercents(coveragePercents)
                .build();
    }
}

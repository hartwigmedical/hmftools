package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.GENOME_TERRITORY_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.MAD_COVERAGE_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.MEAN_COVERAGE_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.MEDIAN_COVERAGE_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.PCT_EXC_ADAPTER_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.PCT_EXC_TOTAL_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.SD_COVERAGE_COLUMN;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;

public final class MetricsWriter
{
    public static void writeResults(final CombinedStats combinedStats, final MetricsConfig config)
    {
        if(config.WriteOldStyle)
        {
            writeOldStyleFile(combinedStats.coverageMetrics(), config);
        }
        else
        {
            writeMetrics(combinedStats.coverageMetrics(), config);
            writeCoverageFrequency(combinedStats.coverageMetrics(), config);
            writeFragmentLengths(combinedStats.fragmentLengths(), config);
        }
    }

    private static final List<Integer> COVERAGE_LEVELS = Lists.newArrayList(1, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100);

    private static String metricsHeader()
    {
        StringJoiner header = new StringJoiner(TSV_DELIM);

        header.add(GENOME_TERRITORY_COLUMN);
        header.add(MEAN_COVERAGE_COLUMN);
        header.add(SD_COVERAGE_COLUMN);
        header.add(MEDIAN_COVERAGE_COLUMN);
        header.add(MAD_COVERAGE_COLUMN);

        header.add(PCT_EXC_ADAPTER_COLUMN);

        for(FilterType type : FilterType.values())
        {
            if(type == FilterType.UNFILTERED)
                continue;

            header.add(type.fileColumn());
        }

        header.add(PCT_EXC_TOTAL_COLUMN);

        for(Integer coverage : COVERAGE_LEVELS)
        {
            header.add(format("PCT_%dX", coverage));
        }

        return header.toString();
    }

    private static String metricsData(final CoverageMetrics metrics)
    {
        StringJoiner tsvData = new StringJoiner(TSV_DELIM);

        final Statistics statistics = metrics.statistics();

        long genomeTerritory = metrics.genomeTerritory();

        tsvData.add(String.valueOf(genomeTerritory));
        tsvData.add(format("%.3f", statistics.Mean));
        tsvData.add(format("%.3f", statistics.StandardDeviation));
        tsvData.add(format("%.0f", statistics.Median));
        tsvData.add(format("%d", statistics.MedianAbsoluteDeviation));

        tsvData.add(format("%.5f", 0.0));

        double filteredPercTotal = 0;
        for(FilterType type : FilterType.values())
        {
            if(type == FilterType.UNFILTERED)
                continue;

            double typePercent = metrics.calcFilteredPercentage(type);
            filteredPercTotal += typePercent;
            tsvData.add(format("%.5f", typePercent));
        }

        tsvData.add(format("%.5f", filteredPercTotal));

        for(Integer coverage : COVERAGE_LEVELS)
        {
            tsvData.add(format("%.5f", metrics.calcCoverageFrequency(coverage, genomeTerritory)));
        }

        // HET SNP values are uncalculated

        return tsvData.toString();
    }

    private static void writeMetrics(final CoverageMetrics metrics, final MetricsConfig config)
    {
        try
        {
            // write summary metrics
            String filename = config.formFilename("metrics");
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write(metricsHeader());
            writer.newLine();
            writer.write(metricsData(metrics));
            writer.newLine();

            writer.close();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write metrics file: {}", e.toString());
        }
    }

    private static void writeCoverageFrequency(final CoverageMetrics metrics, final MetricsConfig config)
    {
        try
        {
            // write coverage frequency for unfiltered aligned bases
            String filename = config.formFilename("coverage");
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Coverage\tCount");
            writer.newLine();

            // collapse for long distributions
            if(metrics.CoverageFrequency.length > 3000)
            {
                int currentCoverage = 0;
                long coverageTotal = metrics.CoverageFrequency[0];

                for(int i = 1; i < metrics.CoverageFrequency.length; ++i)
                {
                    int nextCoverage = CoverageMetrics.getCoverageBucket(i);

                    if(nextCoverage > currentCoverage)
                    {
                        writer.write(format("%d\t%d", currentCoverage, coverageTotal));
                        writer.newLine();

                        currentCoverage = nextCoverage;
                        coverageTotal = 0;
                    }

                    coverageTotal += metrics.CoverageFrequency[i];
                }

                // write the last entry or bucket's data
                writer.write(format("%d\t%d", currentCoverage, coverageTotal));
                writer.newLine();
            }
            else
            {
                for(int i = 0; i < metrics.CoverageFrequency.length; ++i)
                {
                    writer.write(format("%d\t%d", i, metrics.CoverageFrequency[i]));
                    writer.newLine();
                }

            }

            writer.close();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write coverage frequency file: {}", e.toString());
        }
    }

    private static void writeFragmentLengths(final FragmentLengths fragmentLengths, final MetricsConfig config)
    {
        try
        {
            // write summary metrics
            String filename = config.formFilename("frag_lengths");
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("FragmentLength\tCount");
            writer.newLine();

            for(LengthFrequency lengthFrequency : fragmentLengths.lengthFrequencies())
            {
                writer.write(String.format("%d\t%d", lengthFrequency.Length, lengthFrequency.Frequency));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write frag lengths file: {}", e.toString());
        }
    }

    private static void writeOldStyleFile(final CoverageMetrics metrics, final MetricsConfig config)
    {
        try
        {
            // write summary metrics
            String filename;

            if(config.SampleId != null && !config.SampleId.isEmpty())
            {
                filename = WGSMetricsFile.generateFilename(config.OutputDir, config.SampleId);
            }
            else
            {
                String filePrefix = config.BamFile.substring(0, config.BamFile.indexOf(".bam"));
                filename = filePrefix + WGSMetricsFile.FILE_EXTENSION;
            }

            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("## METRICS CLASS");
            writer.newLine();

            writer.write(metricsHeader());
            writer.newLine();
            writer.write(metricsData(metrics));
            writer.newLine();
            writer.newLine();

            writer.write("## HISTOGRAM");
            writer.newLine();

            writer.write("coverage" + TSV_DELIM + "high_quality_coverage_count");
            writer.newLine();

            for(int i = 0; i < metrics.CoverageFrequency.length; ++i)
            {
                writer.write(format("%d\t%d", i, metrics.CoverageFrequency[i]));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write metrics file: {}", e.toString());
        }

    }

}

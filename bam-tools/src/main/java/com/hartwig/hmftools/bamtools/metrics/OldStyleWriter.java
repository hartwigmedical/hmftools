package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.metrics.FlagQCStats.flagStatsPercentages;
import static com.hartwig.hmftools.bamtools.metrics.MetricsWriter.COVERAGE_LEVELS;
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
import java.util.StringJoiner;

import com.hartwig.hmftools.common.metrics.WGSMetricsFile;

public final class OldStyleWriter
{
    public static void writeFlagStats(final FlagStats flagStats, final MetricsConfig config)
    {
        try
        {
            // write flag stats
            String filename = config.OutputDir + config.SampleId + ".flagstat";

            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write(String.format("%s in total (QC-passed reads + QC-failed reads)", flagStats.getTotal().toString()));
            writer.newLine();

            writer.write(String.format("%s primary", flagStats.getPrimary().toString()));
            writer.newLine();

            writer.write(String.format("%s secondary", flagStats.getSecondary().toString()));
            writer.newLine();

            writer.write(String.format("%s supplementary", flagStats.getSupp().toString()));
            writer.newLine();

            writer.write(String.format("%s duplicates", flagStats.getDuplicate().toString()));
            writer.newLine();

            writer.write(String.format("%s primary duplicates", flagStats.getPrimaryDuplicate().toString()));
            writer.newLine();

            writer.write(String.format("%s mapped %s",
                    flagStats.getMapped().toString(), flagStatsPercentages(flagStats.getMapped(), flagStats.getTotal())));
            writer.newLine();

            writer.write(String.format("%s primary mapped %s",
                    flagStats.getPrimaryMapped().toString(), flagStatsPercentages(flagStats.getPrimaryMapped(), flagStats.getPrimary())));
            writer.newLine();

            writer.write(String.format("%s paired in sequencing", flagStats.getPaired().toString()));
            writer.newLine();

            writer.write(String.format("%s read1", flagStats.getRead1().toString()));
            writer.newLine();

            writer.write(String.format("%s read2", flagStats.getRead2().toString()));
            writer.newLine();

            writer.write(String.format(
                    "%s properly paired %s",
                    flagStats.getProperlyPaired().toString(),
                    flagStatsPercentages(flagStats.getProperlyPaired(), flagStats.getPaired())));
            writer.newLine();

            writer.write(String.format("%s with itself and mate mapped", flagStats.getPairMapped().toString()));
            writer.newLine();

            writer.write(String.format("%s singletons %s",
                    flagStats.getSingleton().toString(),
                    flagStatsPercentages(flagStats.getSingleton(), flagStats.getPaired())));
            writer.newLine();

            writer.write(String.format("%s with mate mapped to a different chr", flagStats.getInterChrPairMapped().toString()));
            writer.newLine();

            writer.write(String.format("%s with mate mapped to a different chr (mapQ>=5)", flagStats.getInterChrPairMapQGE5().toString()));
            writer.newLine();

            writer.close();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write flag stats file: {}", e.toString());
        }
    }

    protected static void writeMetricsSummary(final CoverageMetrics metrics, final MetricsConfig config)
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

}

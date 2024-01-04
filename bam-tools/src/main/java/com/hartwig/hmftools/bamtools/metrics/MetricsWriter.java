package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.metrics.FilterType.DUPLICATE;
import static com.hartwig.hmftools.bamtools.metrics.FilterType.LOW_BASE_QUAL;
import static com.hartwig.hmftools.bamtools.metrics.FilterType.LOW_MAP_QUAL;
import static com.hartwig.hmftools.bamtools.metrics.FilterType.MATE_UNMAPPED;
import static com.hartwig.hmftools.bamtools.metrics.FilterType.MAX_COVERAGE;
import static com.hartwig.hmftools.bamtools.metrics.FilterType.OVERLAPPED;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.BamMetricsSummary;
import com.hartwig.hmftools.common.metrics.ImmutableBamMetricsSummary;
import com.hartwig.hmftools.common.metrics.TargetRegionMetrics;

public final class MetricsWriter
{
    public static void writeResults(final CombinedStats combinedStats, final MetricsConfig config)
    {
        if(config.WriteOldStyle)
        {
            OldStyleWriter.writeMetricsSummary(combinedStats.coverageMetrics(), config);
            OldStyleWriter.writeFlagStats(combinedStats.flagStats(), config); // to be redesigned
        }
        else
        {
            writeMetrics(combinedStats.coverageMetrics(), combinedStats.readCounts(), config);
            writeCoverageFrequency(combinedStats.coverageMetrics(), config);
            writeFragmentLengths(combinedStats.fragmentLengths(), config);
            writeTargetRegionStats(combinedStats.targetRegions(), config);

            // TODO - write new format flag stats, or merge into metrics summary
        }
    }

    protected static final List<Integer> COVERAGE_LEVELS = Lists.newArrayList(1, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100);

    private static void writeMetrics(final CoverageMetrics metrics, final ReadCounts readCounts, final MetricsConfig config)
    {
        try
        {
            // write summary metrics
            String filename = BamMetricsSummary.generateFilename(config.OutputDir, config.SampleId);

            List<Integer> coverageLevels = Lists.newArrayList();
            List<Double> coveragePercents = Lists.newArrayList();

            long totalBaseCount = metrics.genomeTerritory();

            for(Integer coverage : COVERAGE_LEVELS)
            {
                coverageLevels.add(coverage);
                coveragePercents.add(metrics.calcCoverageFrequency(coverage, totalBaseCount));
            }

            BamMetricsSummary summary = ImmutableBamMetricsSummary.builder()
                    .totalRegionBases(totalBaseCount)
                    .totalReads(readCounts.TotalReads)
                    .duplicateReads(readCounts.Duplicates)
                    .dualStrandReads(readCounts.DualStrand)
                    .meanCoverage(metrics.statistics().Mean)
                    .sdCoverage(metrics.statistics().StandardDeviation)
                    .medianCoverage((int)round(metrics.statistics().Median))
                    .madCoverage(metrics.statistics().MedianAbsoluteDeviation)
                    .lowMapQualPercent(metrics.calcFilteredPercentage(LOW_MAP_QUAL))
                    .duplicatePercent(metrics.calcFilteredPercentage(DUPLICATE))
                    .unpairedPercent(metrics.calcFilteredPercentage(MATE_UNMAPPED))
                    .lowBaseQualPercent(metrics.calcFilteredPercentage(LOW_BASE_QUAL))
                    .overlappingReadPercent(metrics.calcFilteredPercentage(OVERLAPPED))
                    .cappedCoveragePercent(metrics.calcFilteredPercentage(MAX_COVERAGE))
                    .coverageLevels(coverageLevels)
                    .coveragePercents(coveragePercents)
                    .build();

            summary.write(filename);
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

            // String filename = config.formFilename("coverage");
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
            String filename = config.formFilename("frag_length");

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

    private static void writeTargetRegionStats(final List<TargetRegionStats> targetRegionStats, final MetricsConfig config)
    {
        if(config.TargetRegions.isEmpty())
            return;

        try
        {
            // write summary metrics
            String filename = TargetRegionMetrics.generateFilename(config.OutputDir, config.SampleId);
            TargetRegionMetrics.write(filename, targetRegionStats.stream().map(x -> x.convert()).collect(Collectors.toList()));
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write target regions file: {}", e.toString());
        }
    }
}

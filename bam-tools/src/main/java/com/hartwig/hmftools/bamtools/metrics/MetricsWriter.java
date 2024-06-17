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
import static com.hartwig.hmftools.bamtools.metrics.FlagQCStats.flagStatsPercentages;
import static com.hartwig.hmftools.bamtools.metrics.FlagStatType.INTER_CHR_PAIR_MAPPED;
import static com.hartwig.hmftools.bamtools.metrics.FlagStatType.INTER_CHR_PAIR_MAP_QUAL_GE5;
import static com.hartwig.hmftools.bamtools.metrics.FlagStatType.MAPPED;
import static com.hartwig.hmftools.bamtools.metrics.FlagStatType.PAIRED;
import static com.hartwig.hmftools.bamtools.metrics.FlagStatType.PAIR_MAPPED;
import static com.hartwig.hmftools.bamtools.metrics.FlagStatType.PRIMARY;
import static com.hartwig.hmftools.bamtools.metrics.FlagStatType.PRIMARY_DUPLICATE;
import static com.hartwig.hmftools.bamtools.metrics.FlagStatType.PRIMARY_MAPPED;
import static com.hartwig.hmftools.bamtools.metrics.FlagStatType.PROPERLY_PAIRED;
import static com.hartwig.hmftools.bamtools.metrics.FlagStatType.SECONDARY;
import static com.hartwig.hmftools.bamtools.metrics.FlagStatType.SINGLETON;
import static com.hartwig.hmftools.bamtools.metrics.FlagStatType.SUPPLEMENTARY;
import static com.hartwig.hmftools.bamtools.metrics.FlagStatType.TOTAL;
import static com.hartwig.hmftools.bamtools.metrics.OffTargetFragments.writeOverlapCounts;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.BamMetricsSummary;
import com.hartwig.hmftools.common.metrics.ImmutableBamMetricsSummary;

public class MetricsWriter
{
    private BufferedWriter mTargetRegionsWriter;
    private BufferedWriter mOffTargetHighFragmentOverlapWriter;

    public MetricsWriter(final MetricsConfig config)
    {
        mTargetRegionsWriter = !config.TargetRegions.isEmpty() ? TargetRegionStats.initialiseWriter(config) : null;

        mOffTargetHighFragmentOverlapWriter = !config.TargetRegions.isEmpty() && config.HighFragmentOverlapThreshold > 0 ?
                OffTargetFragments.initialiseEnrichedRegionWriter(config) : null;
    }

    public BufferedWriter targetRegionsWriter() { return mTargetRegionsWriter; }
    public BufferedWriter offTargetHighFragmentOverlapWriter() { return mOffTargetHighFragmentOverlapWriter; }

    public void close()
    {
        closeBufferedWriter(mTargetRegionsWriter);
        closeBufferedWriter(mOffTargetHighFragmentOverlapWriter);
    }

    public static void writeResults(final CombinedStats combinedStats, final MetricsConfig config)
    {
        writeMetrics(combinedStats.coverageMetrics(), combinedStats.readCounts(), config);
        writeCoverageFrequency(combinedStats.coverageMetrics(), config);
        writeFragmentLengths(combinedStats.fragmentLengths(), config);
        writeFlagCounts(combinedStats.flagStats(), config);

        if(config.WriteOffTarget)
            writeOverlapCounts(config, combinedStats.offTargetOverlapCounts());
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
                    .totalReads(readCounts.Total)
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

    private static void writeFlagCounts(final FlagStats flagStats, final MetricsConfig config)
    {
        try
        {
            // write flag stats
            String filename = config.formFilename("flag_counts");

            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write(String.format("%s in total (QC-passed reads + QC-failed reads)", flagStats.statAsString(TOTAL)));
            writer.newLine();

            writer.write(String.format("%s primary", flagStats.statAsString(PRIMARY)));
            writer.newLine();

            writer.write(String.format("%s secondary", flagStats.statAsString(SECONDARY)));
            writer.newLine();

            writer.write(String.format("%s supplementary", flagStats.statAsString(SUPPLEMENTARY)));
            writer.newLine();

            writer.write(String.format("%s duplicates", flagStats.statAsString(FlagStatType.DUPLICATE)));
            writer.newLine();

            writer.write(String.format("%s primary duplicates", flagStats.statAsString(PRIMARY_DUPLICATE)));
            writer.newLine();

            writer.write(String.format("%s mapped %s",
                    flagStats.statAsString(MAPPED), flagStatsPercentages(flagStats.getStat(MAPPED), flagStats.getStat(TOTAL))));
            writer.newLine();

            writer.write(String.format("%s primary mapped %s",
                    flagStats.statAsString(PRIMARY_MAPPED), flagStatsPercentages(flagStats.getStat(PRIMARY_MAPPED), flagStats.getStat(PRIMARY))));
            writer.newLine();

            writer.write(String.format("%s paired in sequencing", flagStats.statAsString(PAIRED)));
            writer.newLine();

            writer.write(String.format("%s read1", flagStats.statAsString(FlagStatType.READ1)));
            writer.newLine();

            writer.write(String.format("%s read2", flagStats.statAsString(FlagStatType.READ2)));
            writer.newLine();

            writer.write(String.format(
                    "%s properly paired %s",
                    flagStats.statAsString(PROPERLY_PAIRED),
                    flagStatsPercentages(flagStats.getStat(PROPERLY_PAIRED), flagStats.getStat(PAIRED))));
            writer.newLine();

            writer.write(String.format("%s with itself and mate mapped", flagStats.statAsString(PAIR_MAPPED)));
            writer.newLine();

            writer.write(String.format("%s singletons %s",
                    flagStats.statAsString(SINGLETON),
                    flagStatsPercentages(flagStats.getStat(SINGLETON), flagStats.getStat(PAIRED))));
            writer.newLine();

            writer.write(String.format("%s with mate mapped to a different chr", flagStats.statAsString(INTER_CHR_PAIR_MAPPED)));
            writer.newLine();

            writer.write(String.format("%s with mate mapped to a different chr (mapQ>=5)", flagStats.statAsString(INTER_CHR_PAIR_MAP_QUAL_GE5)));
            writer.newLine();

            writer.close();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write flag stats file: {}", e.toString());
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
}

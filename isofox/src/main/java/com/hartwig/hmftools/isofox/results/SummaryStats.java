package com.hartwig.hmftools.isofox.results;

import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.isofox.common.FragmentSizeCalcs;

import org.immutables.value.Value;

@Value.Immutable
public abstract class SummaryStats
{
    public abstract String version();

    public abstract int totalFragmentCount();

    public abstract int readLength();

    // 5th, 50th and 95th percentile intronic fragment length
    public abstract double fragmentLength5thPercent();
    public abstract double fragmentLength50thPercent();
    public abstract double fragmentLength95thPercent();

    // proportion of fragments in 7 highly expressed genes
    public abstract double enrichedGenePercent();

    // Median GC (excluding 7 highly expressed genes)
    public abstract double medianGCRatio();

    public static SummaryStats createSummaryStats(
            int totalFragmentCount, int enrichedGeneFragCount, double medianGCRatio, final List<int[]> fragmentLengths, int maxReadLength)
    {
        final VersionInfo version = new VersionInfo("isofox.version");
        double enrichedGenePercent = totalFragmentCount > 0 ? enrichedGeneFragCount / (double)totalFragmentCount : 0;

        final List<Double> fragLengths = FragmentSizeCalcs.calcPercentileData(fragmentLengths, Lists.newArrayList(0.05, 0.5, 0.95));

        return ImmutableSummaryStats.builder()
                .version(version.version())
                .totalFragmentCount(totalFragmentCount)
                .readLength(maxReadLength)
                .enrichedGenePercent(enrichedGenePercent)
                .medianGCRatio(medianGCRatio)
                .fragmentLength5thPercent(fragLengths.get(0))
                .fragmentLength50thPercent(fragLengths.get(1))
                .fragmentLength95thPercent(fragLengths.get(2))
                .build();
    }

    public static String csvHeader()
    {
        return "SampleId,Version,TotalFragments,ReadLength,FragLength5th,FragLength50th,FragLength95th"
                + ",EnrichedGenePercent,MedianGCRatio";
    }

    public String toCsv(final String sampleId)
    {
        return new StringJoiner(DELIMITER)
                .add(sampleId)
                .add(version())
                .add(String.valueOf(totalFragmentCount()))
                .add(String.valueOf(readLength()))
                .add(String.format("%.0f", fragmentLength5thPercent()))
                .add(String.format("%.0f", fragmentLength50thPercent()))
                .add(String.format("%.0f", fragmentLength95thPercent()))
                .add(String.format("%.3f", enrichedGenePercent()))
                .add(String.format("%.3f", medianGCRatio()))
                .toString();
    }

    public static SummaryStats fromCsv(final String input)
    {
        final String items[] = input.split(",");
        int index = 1;

        return ImmutableSummaryStats.builder()
                .version(items[index++])
                .totalFragmentCount(Integer.parseInt(items[index++]))
                .readLength(Integer.parseInt(items[index++]))
                .fragmentLength5thPercent(Double.parseDouble(items[index++]))
                .fragmentLength50thPercent(Double.parseDouble(items[index++]))
                .fragmentLength95thPercent(Double.parseDouble(items[index++]))
                .enrichedGenePercent(Double.parseDouble(items[index++]))
                .medianGCRatio(Double.parseDouble(items[index++]))
                .build();
    }

    public static SummaryStats loadFile(final Path filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(filename);
            if(lines.size() != 2)
            {
                ISF_LOGGER.error("failed to load file({}: invalid line count({})", filename.toString(), lines.size());
                return null;
            }

            return fromCsv(lines.get(1));
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load file({}): {}", filename.toString(), e.toString());
            return null;
        }
    }

}

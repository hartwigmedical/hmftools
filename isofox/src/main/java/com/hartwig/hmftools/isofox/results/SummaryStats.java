package com.hartwig.hmftools.isofox.results;

import java.util.List;

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

}

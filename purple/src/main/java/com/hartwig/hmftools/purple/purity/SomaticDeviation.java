package com.hartwig.hmftools.purple.purity;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.cache.RemovalListener;
import com.google.common.cache.RemovalNotification;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public enum SomaticDeviation implements RemovalListener<Double, Integer>
{
    INSTANCE;

    private static final int TRIALS = 10_000;
    private static final Logger LOGGER = LogManager.getLogger(SomaticDeviation.class);
    private final LoadingCache<Double, Integer> maxConceivableCache;
    private boolean errorLogged = false;

    SomaticDeviation()
    {
        this.maxConceivableCache =
                CacheBuilder.newBuilder().maximumSize(50000).removalListener(this).build(new CacheLoader<Double, Integer>()
                {

                    @Override
                    public Integer load(@NotNull final Double p)
                    {
                        final BinomialDistribution dist = new BinomialDistribution(TRIALS, Math.min(1, p));
                        return dist.inverseCumulativeProbability(0.999);
                    }
                });
    }

    public double deviationFromMax(@NotNull final PurityAdjuster purityAdjuster, @NotNull final FittedRegion region,
            @NotNull final SomaticVariant variant)
    {
        double normalCopyNumber = purityAdjuster.germlineCopyNumber(region.chromosome());
        double constrainedMajorAllelePloidy = Math.max(0, region.majorAlleleCopyNumber());
        double constrainedTumorCopyNumber = Math.max(0, region.tumorCopyNumber());

        return deviationFromMax(purityAdjuster, normalCopyNumber, variant, constrainedTumorCopyNumber, constrainedMajorAllelePloidy);
    }

    @VisibleForTesting
    double deviationFromMax(@NotNull final PurityAdjuster purityAdjuster, double normalCopyNumber, @NotNull final AllelicDepth depth,
            double tumorCopyNumber, double tumorMajorAllelePloidy)
    {
        double maxConceivablePloidy =
                maxConceivablePloidy(purityAdjuster, normalCopyNumber, depth, tumorCopyNumber, tumorMajorAllelePloidy);
        double somaticPloidy = purityAdjuster.purityAdjustedPloidy(normalCopyNumber, 0, tumorCopyNumber, depth.alleleFrequency());

        return Math.max(0, somaticPloidy - maxConceivablePloidy);
    }

    @VisibleForTesting
    double maxConceivablePloidy(@NotNull final PurityAdjuster purityAdjuster, double normalCopyNumber, @NotNull final AllelicDepth depth,
            double tumorCopyNumber, double tumorMajorAllelePloidy)
    {
        final int maxConceivableReads =
                maxConceivableReads(purityAdjuster, normalCopyNumber, depth, tumorCopyNumber, tumorMajorAllelePloidy);
        final double maxConceivableVAF = 1d * maxConceivableReads / depth.totalReadCount();

        return purityAdjuster.purityAdjustedPloidy(normalCopyNumber, 0, tumorCopyNumber, maxConceivableVAF);
    }

    @VisibleForTesting
    int maxConceivableReads(@NotNull final PurityAdjuster purityAdjuster, double normalCopyNumber, @NotNull final AllelicDepth depth,
            double tumorCopyNumber, double tumorMajorAllelePloidy)
    {
        double expectedVAF = purityAdjuster.expectedFrequency(normalCopyNumber, 0, tumorCopyNumber, tumorMajorAllelePloidy);
        double p = 1d * Math.round(expectedVAF * depth.totalReadCount() * 100) / 100 / TRIALS;
        return maxConceivableCache.getUnchecked(p);
    }

    @Override
    public void onRemoval(@NotNull final RemovalNotification<Double, Integer> removalNotification)
    {
        if(!errorLogged)
        {
            errorLogged = true;
            LOGGER.warn("Somatic deviation cache limit exceeded. This indicates a potential performance issue but will otherwise not effect any calculations.");
        }
    }
}

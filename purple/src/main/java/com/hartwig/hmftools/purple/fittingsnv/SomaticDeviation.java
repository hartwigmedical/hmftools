package com.hartwig.hmftools.purple.fittingsnv;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.cache.RemovalListener;
import com.google.common.cache.RemovalNotification;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

import org.apache.commons.math3.distribution.BinomialDistribution;

public enum SomaticDeviation implements RemovalListener<Double, Integer>
{
    INSTANCE;

    private static final int TRIALS = 10_000;
    private final LoadingCache<Double, Integer> mMaxConceivableCache;
    private boolean mErrorLogged = false;

    SomaticDeviation()
    {
        mMaxConceivableCache =
                CacheBuilder.newBuilder().maximumSize(50000).removalListener(this).build(new CacheLoader<Double, Integer>()
                {

                    @Override
                    public Integer load(final Double p)
                    {
                        final BinomialDistribution dist = new BinomialDistribution(TRIALS, Math.min(1, p));
                        return dist.inverseCumulativeProbability(0.999);
                    }
                });
    }

    public double deviationFromMax(
            final PurityAdjuster purityAdjuster, final String chromosome, double majorAlleleCopyNumber, double tumorCopyNumber,
            final SomaticVariant variant)
    {
        double normalCopyNumber = purityAdjuster.germlineCopyNumber(chromosome);
        double constrainedMajorAllelePloidy = Math.max(0, majorAlleleCopyNumber);
        double constrainedTumorCopyNumber = Math.max(0, tumorCopyNumber);

        return deviationFromMax(
                purityAdjuster, normalCopyNumber, variant.tumorAlleleDepth(), constrainedTumorCopyNumber, constrainedMajorAllelePloidy);
    }

    @VisibleForTesting
    public double deviationFromMax(
            final PurityAdjuster purityAdjuster, double normalCopyNumber, final AllelicDepth depth,
            double tumorCopyNumber, double tumorMajorAllelePloidy)
    {
        double maxConceivablePloidy =
                maxConceivablePloidy(purityAdjuster, normalCopyNumber, depth, tumorCopyNumber, tumorMajorAllelePloidy);
        double somaticPloidy = purityAdjuster.purityAdjustedPloidy(normalCopyNumber, 0, tumorCopyNumber, depth.alleleFrequency());

        return Math.max(0, somaticPloidy - maxConceivablePloidy);
    }

    @VisibleForTesting
    public double maxConceivablePloidy(
            final PurityAdjuster purityAdjuster, double normalCopyNumber, final AllelicDepth depth,
            double tumorCopyNumber, double tumorMajorAllelePloidy)
    {
        final int maxConceivableReads =
                maxConceivableReads(purityAdjuster, normalCopyNumber, depth, tumorCopyNumber, tumorMajorAllelePloidy);
        final double maxConceivableVAF = 1d * maxConceivableReads / depth.TotalReadCount;

        return purityAdjuster.purityAdjustedPloidy(normalCopyNumber, 0, tumorCopyNumber, maxConceivableVAF);
    }

    @VisibleForTesting
    public int maxConceivableReads(
            final PurityAdjuster purityAdjuster, double normalCopyNumber, final AllelicDepth depth,
            double tumorCopyNumber, double tumorMajorAllelePloidy)
    {
        double expectedVAF = purityAdjuster.expectedFrequency(normalCopyNumber, 0, tumorCopyNumber, tumorMajorAllelePloidy);
        double p = 1d * Math.round(expectedVAF * depth.TotalReadCount * 100) / 100 / TRIALS;
        return mMaxConceivableCache.getUnchecked(p);
    }

    @Override
    public void onRemoval(final RemovalNotification<Double, Integer> removalNotification)
    {
        if(!mErrorLogged)
        {
            mErrorLogged = true;
            PPL_LOGGER.debug("somatic deviation cache limit exceeded");
        }
    }
}

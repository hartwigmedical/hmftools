package com.hartwig.hmftools.common.purple.purity;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.jetbrains.annotations.NotNull;

public class SomaticDeviation {

    private static final int TRIALS = 100_000;
    private final PurityAdjuster purityAdjuster;

    public SomaticDeviation(final PurityAdjuster purityAdjuster) {
        this.purityAdjuster = purityAdjuster;
    }

    public double deviationFromMax(@NotNull final FittedRegion region, @NotNull final SomaticVariant variant) {
        int normalCopyNumber = purityAdjuster.typicalCopyNumber(region.chromosome());
        double constrainedMajorAllelePloidy = Math.max(0, region.majorAllelePloidy());
        double constrainedTumorCopyNumber = Math.max(0, region.tumorCopyNumber());

        return deviationFromMax(normalCopyNumber, variant, constrainedTumorCopyNumber, constrainedMajorAllelePloidy);
    }

    @VisibleForTesting
    double deviationFromMax(int normalCopyNumber, @NotNull final AllelicDepth depth, double tumorCopyNumber,
            double tumorMajorAllelePloidy) {
        double maxConceivablePloidy = maxConceivablePloidy(normalCopyNumber, depth, tumorCopyNumber, tumorMajorAllelePloidy);
        double somaticPloidy = purityAdjuster.purityAdjustedPloidy(normalCopyNumber, 0, tumorCopyNumber, depth.alleleFrequency());

        return Math.max(0, somaticPloidy - maxConceivablePloidy);
    }

    @VisibleForTesting
    double maxConceivablePloidy(int normalCopyNumber, @NotNull final AllelicDepth depth, double tumorCopyNumber,
            double tumorMajorAllelePloidy) {
        final int maxConceivableReads = maxConceivableReads(normalCopyNumber, depth, tumorCopyNumber, tumorMajorAllelePloidy);
        final double maxConceivableVAF = 1d * maxConceivableReads / depth.totalReadCount();
        final double maxConceivablePloidy = purityAdjuster.purityAdjustedPloidy(normalCopyNumber, 0, tumorCopyNumber, maxConceivableVAF);

        return maxConceivablePloidy;
    }

    @VisibleForTesting
    int maxConceivableReads(int normalCopyNumber, @NotNull final AllelicDepth depth, double tumorCopyNumber,
            double tumorMajorAllelePloidy) {
        double expectedVAF = purityAdjuster.expectedFrequency(normalCopyNumber, 0, tumorCopyNumber, tumorMajorAllelePloidy);

        final BinomialDistribution dist = new BinomialDistribution(TRIALS, Math.min(1, depth.totalReadCount() * expectedVAF / TRIALS));
        final int maxConceivableReads = dist.inverseCumulativeProbability(0.99);
        return maxConceivableReads;
    }

}

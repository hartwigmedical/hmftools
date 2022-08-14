package com.hartwig.hmftools.purple.purity;

import java.util.Collection;
import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

public final class SomaticPenaltyFactory
{
    public static double calcPenalty(
            final PurityAdjuster purityAdjuster, final Collection<ObservedRegion> regions, final List<SomaticVariant> variants)
    {
        final SomaticDeviation somaticDeviation = SomaticDeviation.INSTANCE;

        final GenomePositionSelector<SomaticVariant> variantSelector = GenomePositionSelectorFactory.create(variants);
        double score = 0;
        int variantCount = 0;

        for(ObservedRegion region : regions)
        {
            SomaticVariantConsumer consumer = new SomaticVariantConsumer(purityAdjuster, somaticDeviation, region);
            variantSelector.select(region, consumer);
            score += consumer.score();
            variantCount += consumer.variants();
        }

        return variantCount == 0 ? 0 : score / variantCount;
    }

    private static class SomaticVariantConsumer implements Consumer<SomaticVariant>
    {
        final PurityAdjuster mPurityAdjuster;
        private final SomaticDeviation mSomaticDeviation;
        private final ObservedRegion mRegion;
        private double mScore;
        private int mVariants;

        private SomaticVariantConsumer(
                final PurityAdjuster purityAdjuster, final SomaticDeviation somaticDeviation, final ObservedRegion region)
        {
            mPurityAdjuster = purityAdjuster;
            mSomaticDeviation = somaticDeviation;
            mRegion = region;
        }

        @Override
        public void accept(final SomaticVariant variant)
        {
            mScore += mSomaticDeviation.deviationFromMax(mPurityAdjuster, mRegion, variant);
            mVariants++;
        }

        public int variants()
        {
            return mVariants;
        }

        public double score()
        {
            return mScore;
        }
    }
}

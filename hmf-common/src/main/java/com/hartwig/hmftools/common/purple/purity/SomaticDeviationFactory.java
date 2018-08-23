package com.hartwig.hmftools.common.purple.purity;

import java.util.Collection;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class SomaticDeviationFactory {

    private final SomaticDeviation somaticDeviation;

    protected SomaticDeviationFactory(final PurityAdjuster purityAdjuster) {
        somaticDeviation = new SomaticDeviation(purityAdjuster);
    }

    public double deviation(@NotNull Collection<FittedRegion> regions, Collection<SomaticVariant> variants) {
        final GenomePositionSelector<SomaticVariant> variantSelector = GenomePositionSelectorFactory.create(variants);
        double score = 0;

        for (FittedRegion region : regions) {
            SomaticVariantConsumer consumer = new SomaticVariantConsumer(region);
            variantSelector.select(region, consumer);
            score += consumer.score();
        }

        return score / variants.size();
    }

    private class SomaticVariantConsumer implements Consumer<SomaticVariant> {

        private final FittedRegion region;
        private double score;

        private SomaticVariantConsumer(final FittedRegion region) {
            this.region = region;
        }

        @Override
        public void accept(final SomaticVariant variant) {
            score += somaticDeviation.deviationFromMax(region, variant);
        }

        public double score() {
            return score;
        }
    }
}

package com.hartwig.hmftools.common.purple.region;

import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.segment.PurpleSegment;

import org.jetbrains.annotations.NotNull;

public class ObservedRegionFactory {

    private final Gender gender;
    private final ObservedRegionStatusFactory statusFactory;

    public ObservedRegionFactory(final Gender gender) {
        this.gender = gender;
        statusFactory = new ObservedRegionStatusFactory(gender);
    }

    @NotNull
    public List<ObservedRegion> combine(@NotNull final List<PurpleSegment> regions, @NotNull final Multimap<String, AmberBAF> bafs,
            @NotNull final Multimap<String, CobaltRatio> ratios) {
        final List<ObservedRegion> result = Lists.newArrayList();

        final GenomePositionSelector<CobaltRatio> cobaltSelector = GenomePositionSelectorFactory.create(ratios);
        final GenomePositionSelector<AmberBAF> bafSelector = GenomePositionSelectorFactory.create(bafs);

        for (final PurpleSegment region : regions) {
            final BAFAccumulator baf = new BAFAccumulator();
            final CobaltAccumulator cobalt = new CobaltAccumulator();

            bafSelector.select(region, baf);
            cobaltSelector.select(region, cobalt);

            double tumorRatio = cobalt.tumorMeanRatio();
            double normalRatio = cobalt.referenceMeanRatio();
            final EnrichedRegion copyNumber = ImmutableEnrichedRegion.builder()
                    .from(region)
                    .bafCount(baf.count())
                    .observedBAF(baf.medianBaf())
                    .observedTumorRatio(tumorRatio)
                    .observedNormalRatio(normalRatio)
                    .ratioSupport(region.ratioSupport())
                    .support(region.support())
                    .observedTumorRatioCount(cobalt.tumorCount())
                    .status(statusFactory.status(region, normalRatio))
                    .build();

            result.add(copyNumber);
        }

        return result;
    }

    private class BAFAccumulator implements Consumer<AmberBAF> {
        private int count;
        final private List<Double> bafs = Lists.newArrayList();

        @Override
        public void accept(final AmberBAF baf) {
            if (HumanChromosome.valueOf(baf).isHomologous(gender) && !Double.isNaN(baf.tumorModifiedBAF())) {
                count++;
                bafs.add(baf.tumorModifiedBAF());
            }
        }

        private int count() {
            return count;
        }

        private double medianBaf() {
            if (count > 0) {
                Collections.sort(bafs);
                return bafs.size() % 2 == 0 ? (bafs.get(count / 2) + bafs.get(count / 2 - 1)) / 2 : bafs.get(count / 2);
            }
            return 0;
        }
    }

    private class CobaltAccumulator implements Consumer<CobaltRatio> {
        private final RatioAccumulator referenceAccumulator = new RatioAccumulator();
        private final RatioAccumulator tumorAccumulator = new RatioAccumulator();

        private double referenceMeanRatio() {
            return referenceAccumulator.meanRatio();
        }

        private double tumorMeanRatio() {
            return tumorAccumulator.meanRatio();
        }

        private int tumorCount() {
            return tumorAccumulator.count();
        }

        @Override
        public void accept(final CobaltRatio ratio) {
            referenceAccumulator.accept(ratio.referenceGCDiploidRatio());
            tumorAccumulator.accept(ratio.tumorGCRatio());
        }
    }

    private class RatioAccumulator implements Consumer<Double> {
        private double sumRatio;
        private int count;

        private double meanRatio() {
            return count > 0 ? sumRatio / count : 0;
        }

        private int count() {
            return count;
        }

        @Override
        public void accept(final Double ratio) {
            if (Doubles.greaterThan(ratio, -1)) {
                count++;
                sumRatio += ratio;
            }
        }
    }
}

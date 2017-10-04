package com.hartwig.hmftools.common.purple.region;

import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.baf.TumorBAF;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.copynumber.freec.FreecStatus;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.segment.PurpleSegment;

import org.jetbrains.annotations.NotNull;

public class ObservedRegionFactory {

    private final Gender gender;

    public ObservedRegionFactory(final Gender gender) {
        this.gender = gender;
    }

    @NotNull
    public List<ObservedRegion> combine(@NotNull final List<PurpleSegment> regions, @NotNull final Multimap<String, TumorBAF> bafs,
            @NotNull final Multimap<String, ReadRatio> tumorRatios, @NotNull final Multimap<String, ReadRatio> normalRatios) {
        final List<ObservedRegion> result = Lists.newArrayList();

        final GenomePositionSelector<ReadRatio> tumorRatioSelector = GenomePositionSelectorFactory.create(tumorRatios);
        final GenomePositionSelector<ReadRatio> normalRatioSelector = GenomePositionSelectorFactory.create(normalRatios);
        final GenomePositionSelector<TumorBAF> bafSelector = GenomePositionSelectorFactory.create(bafs);

        for (final PurpleSegment region : regions) {
            final BAFAccumulator baf = new BAFAccumulator();
            final RatioAccumulator tumorRatio = new RatioAccumulator();
            final RatioAccumulator normalRatio = new RatioAccumulator();

            bafSelector.select(region, baf);
            tumorRatioSelector.select(region, tumorRatio);
            normalRatioSelector.select(region, normalRatio);

            double myTumorRatio = tumorRatio.meanRatio();
            double myNormalRatio = normalRatio.meanRatio();
            final EnrichedRegion copyNumber = ImmutableEnrichedRegion.builder()
                    .from(region)
                    .bafCount(baf.count())
                    .observedBAF(baf.medianBaf())
                    .observedTumorRatio(myTumorRatio)
                    .observedNormalRatio(myNormalRatio)
                    .ratioSupport(region.ratioSupport())
                    .structuralVariantSupport(region.structuralVariantSupport())
                    .observedTumorRatioCount(tumorRatio.count())
                    .status(FreecStatus.fromNormalRatio(gender, region.chromosome(), myNormalRatio))
                    .build();

            result.add(copyNumber);
        }

        return result;
    }

    private class BAFAccumulator implements Consumer<TumorBAF> {
        private int count;
        final private List<Double> bafs = Lists.newArrayList();

        @Override
        public void accept(final TumorBAF baf) {
            if (HumanChromosome.valueOf(baf).isHomologous(gender) && !Double.isNaN(baf.mBaf())) {
                count++;
                bafs.add(baf.mBaf());
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

    private class RatioAccumulator implements Consumer<ReadRatio> {
        private double sumRatio;
        private int count;

        private double meanRatio() {
            return count > 0 ? sumRatio / count : 0;
        }

        private int count() {
            return count;
        }

        @Override
        public void accept(final ReadRatio ratio) {
            if (Doubles.greaterThan(ratio.ratio(), -1)) {
                count++;
                sumRatio += ratio.ratio();
            }
        }
    }
}

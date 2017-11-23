package com.hartwig.hmftools.common.purple.region;

import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.gc.GCProfile;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.segment.PurpleSegment;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.window.Window;

import org.jetbrains.annotations.NotNull;

public class ObservedRegionFactory {

    private final int windowSize;
    private final Gender gender;
    private final ObservedRegionStatusFactory statusFactory;

    public ObservedRegionFactory(final int windowSize, final Gender gender) {
        this.windowSize = windowSize;
        this.gender = gender;
        statusFactory = new ObservedRegionStatusFactory(gender);
    }

    @NotNull
    public List<ObservedRegion> combine(@NotNull final List<PurpleSegment> regions, @NotNull final Multimap<String, AmberBAF> bafs,
            @NotNull final Multimap<String, CobaltRatio> ratios, @NotNull final Multimap<String, GCProfile> gcProfiles) {
        final List<ObservedRegion> result = Lists.newArrayList();

        final GenomePositionSelector<CobaltRatio> cobaltSelector = GenomePositionSelectorFactory.create(ratios);
        final GenomePositionSelector<AmberBAF> bafSelector = GenomePositionSelectorFactory.create(bafs);
        final GenomeRegionSelector<GCProfile> gcSelector = GenomeRegionSelectorFactory.create(gcProfiles);

        for (final PurpleSegment region : regions) {
            final BAFAccumulator baf = new BAFAccumulator();
            final CobaltAccumulator cobalt = new CobaltAccumulator(windowSize, region);
            final GCAccumulator gc = new GCAccumulator(region);

            bafSelector.select(region, baf);
            cobaltSelector.select(region, cobalt);
            gcSelector.select(region, gc);

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
                    .gcContent(gc.averageGCContent())
                    .status(statusFactory.status(region, normalRatio, tumorRatio))
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

    @VisibleForTesting
    static class CobaltAccumulator implements Consumer<CobaltRatio> {

        private final Window window;
        private final GenomeRegion region;

        private final RatioAccumulator referenceAccumulator = new RatioAccumulator();
        private final RatioAccumulator tumorAccumulator = new RatioAccumulator();

        CobaltAccumulator(final int windowSize, final GenomeRegion region) {
            this.window = new Window(windowSize);
            this.region = region;
        }

        double referenceMeanRatio() {
            return referenceAccumulator.meanRatio();
        }

        double tumorMeanRatio() {
            return tumorAccumulator.meanRatio();
        }

        int tumorCount() {
            return tumorAccumulator.count();
        }

        @Override
        public void accept(final CobaltRatio ratio) {
            if (window.end(ratio.position()) <= region.end()) {
                referenceAccumulator.accept(ratio.referenceGCDiploidRatio());
                tumorAccumulator.accept(ratio.tumorGCRatio());
            }
        }
    }

    static private class RatioAccumulator implements Consumer<Double> {
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

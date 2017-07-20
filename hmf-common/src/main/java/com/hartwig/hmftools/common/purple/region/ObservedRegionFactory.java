package com.hartwig.hmftools.common.purple.region;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.copynumber.freec.FreecRatio;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.purple.ratio.GCContent;
import com.hartwig.hmftools.common.purple.segment.PurpleSegment;
import com.hartwig.hmftools.common.variant.GermlineSampleData;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public class ObservedRegionFactory {

    private static final Set<String> HETEROZYGOUS_GENO_TYPES = Sets.newHashSet("0/1", "0|1");

    private final double minRefAlleleFrequency;
    private final double maxRefAlleleFrequency;
    private final long minCombinedDepth;
    private final long maxCombinedDepth;

    public ObservedRegionFactory(final double minRefAlleleFrequency, final double maxRefAlleleFrequency, final long minCombinedDepth,
            final long maxCombinedDepth) {
        this.minRefAlleleFrequency = minRefAlleleFrequency;
        this.maxRefAlleleFrequency = maxRefAlleleFrequency;
        this.minCombinedDepth = minCombinedDepth;
        this.maxCombinedDepth = maxCombinedDepth;
    }

    @NotNull
    public List<ObservedRegion> combine(@NotNull final List<PurpleSegment> regions, @NotNull final List<GermlineVariant> variants,
            @NotNull final List<FreecRatio> tumorRatios, @NotNull final List<FreecRatio> normalRatios,
            @NotNull final Multimap<String, GCContent> gcContents) {
        final List<ObservedRegion> result = Lists.newArrayList();

        final GenomePositionSelector<FreecRatio> tumorRatioSelector = GenomePositionSelectorFactory.create(tumorRatios);
        final GenomePositionSelector<FreecRatio> normalRatioSelector = GenomePositionSelectorFactory.create(normalRatios);
        final GenomePositionSelector<GermlineVariant> variantSelector = GenomePositionSelectorFactory.create(variants);
        final GenomePositionSelector<GCContent> gcContentSelector = GenomePositionSelectorFactory.create(gcContents);

        for (final PurpleSegment region : regions) {
            final BAFAccumulator baf = new BAFAccumulator();
            final RatioAccumulator tumorRatio = new RatioAccumulator();
            final RatioAccumulator normalRatio = new RatioAccumulator();
            final GCContentAccumulator gcContent = new GCContentAccumulator();

            variantSelector.select(region, baf);
            tumorRatioSelector.select(region, tumorRatio);
            normalRatioSelector.select(region, normalRatio);
            gcContentSelector.select(region, gcContent);

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
                    .observedGCContent(gcContent.getAverageGCContent())
                    .observedNonNPercentage(gcContent.getAverageNonNPercentage())
                    .observedMappablePercentage(gcContent.getAverageMappablePercentage())
                    .build();

            result.add(copyNumber);
        }

        return result;
    }

    private class BAFAccumulator implements Consumer<GermlineVariant> {
        private int count;
        final private List<Double> bafs = Lists.newArrayList();

        @Override
        public void accept(final GermlineVariant variant) {
            final GermlineSampleData tumorData = variant.tumorData();
            if (tumorData == null || !HETEROZYGOUS_GENO_TYPES.contains(variant.refData().genoType()) || variant.type() != VariantType.SNP
                    || variant.refData().alleleFrequency() <= minRefAlleleFrequency
                    || variant.refData().alleleFrequency() >= maxRefAlleleFrequency || variant.refData().combinedDepth() <= minCombinedDepth
                    || variant.refData().combinedDepth() >= maxCombinedDepth) {
                return;
            }

            double standardBAF = tumorData.alleleFrequency();
            double modifiedBAF = 0.5 + Math.abs(standardBAF - 0.5);

            count++;
            bafs.add(modifiedBAF);
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

    private class RatioAccumulator implements Consumer<FreecRatio> {
        private double sumRatio;
        private int count;

        private double meanRatio() {
            return count > 0 ? sumRatio / count : 0;
        }

        @Override
        public void accept(final FreecRatio ratio) {
            if (Doubles.greaterThan(ratio.ratio(), -1)) {
                count++;
                sumRatio += ratio.ratio();
            }
        }
    }

    @VisibleForTesting
    static class GCContentAccumulator implements Consumer<GCContent> {

        private int count;
        private double sumGCContent;
        private double sumNonNPercentage;
        private double sumMappablePercentage;

        @Override
        public void accept(final GCContent GCContent) {
            count++;
            if (Doubles.positiveOrZero(GCContent.gcContent())) {
                sumGCContent += GCContent.gcContent() * GCContent.nonNPercentage();
                sumNonNPercentage += GCContent.nonNPercentage();
                sumMappablePercentage += GCContent.mappablePercentage();
            }
        }

        double getAverageGCContent() {
            return Doubles.isZero(sumNonNPercentage) ? 0 : sumGCContent / sumNonNPercentage;
        }

        double getAverageNonNPercentage() {
            return count == 0 ? 0 : sumNonNPercentage / count;
        }

        double getAverageMappablePercentage() {
            return count == 0 ? 0 : sumMappablePercentage / count;
        }
    }
}

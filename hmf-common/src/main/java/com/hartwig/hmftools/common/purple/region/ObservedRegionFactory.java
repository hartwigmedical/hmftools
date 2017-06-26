package com.hartwig.hmftools.common.purple.region;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.copynumber.freec.FreecRatio;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.segment.PurpleSegment;
import com.hartwig.hmftools.common.variant.GermlineSampleData;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.zipper.GenomeZipper;
import com.hartwig.hmftools.common.zipper.GenomeZipperRegionHandler;

import org.jetbrains.annotations.NotNull;

public class ObservedRegionFactory implements GenomeZipperRegionHandler<PurpleSegment> {

    private static final Set<String> HETEROZYGOUS_GENO_TYPES = Sets.newHashSet("0/1", "0|1");

    private final double minRefAlleleFrequency;
    private final double maxRefAlleleFrequency;
    private final long minCombinedDepth;
    private final long maxCombinedDepth;

    private final List<ObservedRegion> result = Lists.newArrayList();
    private final BAFAccumulator baf = new BAFAccumulator();
    private final RatioAccumulator tumorRatio = new RatioAccumulator();
    private final RatioAccumulator normalRatio = new RatioAccumulator();

    public ObservedRegionFactory(final double minRefAlleleFrequency, final double maxRefAlleleFrequency, final long minCombinedDepth,
            final long maxCombinedDepth) {
        this.minRefAlleleFrequency = minRefAlleleFrequency;
        this.maxRefAlleleFrequency = maxRefAlleleFrequency;
        this.minCombinedDepth = minCombinedDepth;
        this.maxCombinedDepth = maxCombinedDepth;
    }

    @NotNull
    public List<ObservedRegion> combine(@NotNull final List<PurpleSegment> regions, @NotNull final List<GermlineVariant> variants,
            @NotNull final List<FreecRatio> tumorRatios, @NotNull final List<FreecRatio> normalRatios) {
        baf.reset();
        tumorRatio.reset();
        normalRatio.reset();
        result.clear();

        final GenomeZipper<PurpleSegment> zipper = new GenomeZipper<>(false, regions, this);
        zipper.addPositions(variants, this::variant);
        zipper.addPositions(tumorRatios, tumorRatio::accumulate);
        zipper.addPositions(normalRatios, normalRatio::accumulate);
        zipper.zip();

        return result;
    }

    @Override
    public void chromosome(@NotNull final String chromosome) {

    }

    @Override
    public void enter(@NotNull final PurpleSegment region) {
        baf.reset();
        tumorRatio.reset();
        normalRatio.reset();
    }

    @Override
    public void exit(@NotNull final PurpleSegment region) {
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
                .build();

        result.add(copyNumber);
    }

    private void variant(@NotNull final GermlineVariant variant) {
        final GermlineSampleData tumorData = variant.tumorData();
        if (tumorData == null || !HETEROZYGOUS_GENO_TYPES.contains(variant.refData().genoType()) || variant.type() != VariantType.SNP
                || variant.refData().alleleFrequency() <= minRefAlleleFrequency
                || variant.refData().alleleFrequency() >= maxRefAlleleFrequency || variant.refData().combinedDepth() <= minCombinedDepth
                || variant.refData().combinedDepth() >= maxCombinedDepth) {
            return;
        }

        double standardBAF = tumorData.alleleFrequency();
        double modifiedBAF = 0.5 + Math.abs(standardBAF - 0.5);
        baf.accumulate(modifiedBAF);
    }

    private class BAFAccumulator {
        private int count;
        final private List<Double> bafs = Lists.newArrayList();

        private void reset() {
            count = 0;
            bafs.clear();
        }

        private void accumulate(final double baf) {
            count++;
            bafs.add(baf);
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

    private class RatioAccumulator {
        private double sumRatio;
        private int count;

        private void accumulate(@NotNull final FreecRatio ratio) {
            if (Doubles.greaterThan(ratio.ratio(), -1)) {
                count++;
                sumRatio += ratio.ratio();
            }
        }

        private void reset() {
            count = 0;
            sumRatio = 0;
        }

        private double meanRatio() {
            return count > 0 ? sumRatio / count : 0;
        }
    }
}

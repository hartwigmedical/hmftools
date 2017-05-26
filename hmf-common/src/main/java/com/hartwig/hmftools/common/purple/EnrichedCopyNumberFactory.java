package com.hartwig.hmftools.common.purple;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.copynumber.freec.FreecRatio;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.GermlineSampleData;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.zipper.GenomeZipper;
import com.hartwig.hmftools.common.zipper.GenomeZipperRegionHandler;

public class EnrichedCopyNumberFactory implements GenomeZipperRegionHandler<GenomeRegion> {

    private static final Set<String> GENO_TYPE = Sets.newHashSet("0/1", "0|1");
    private final double minRefAlleleFrequency;
    private final double maxRefAlleleFrequency;
    private final long minCombinedDepth;
    private final long maxCombinedDepth;

    private final List<EnrichedCopyNumber> result = Lists.newArrayList();
    private final BAFAccumulator baf = new BAFAccumulator();
    private final RatioAccumulator tumorRatio = new RatioAccumulator();
    private final RatioAccumulator normalRatio = new RatioAccumulator();

    public EnrichedCopyNumberFactory(double minRefAlleleFrequency, double maxRefAlleleFrequency, long minCombinedDepth,
            long maxCombinedDepth) {
        this.minRefAlleleFrequency = minRefAlleleFrequency;
        this.maxRefAlleleFrequency = maxRefAlleleFrequency;
        this.minCombinedDepth = minCombinedDepth;
        this.maxCombinedDepth = maxCombinedDepth;
    }

    public List<EnrichedCopyNumber> enrich(List<GenomeRegion> copyNumbers, List<GermlineVariant> variants,
            List<FreecRatio> tumorRatios, List<FreecRatio> normalRatios) {
        baf.reset();
        tumorRatio.reset();
        normalRatio.reset();
        result.clear();

        GenomeZipper<GenomeRegion> zipper = new GenomeZipper<>(copyNumbers, this);
        zipper.addPositions(variants, this::variant);
        zipper.addPositions(tumorRatios, tumorRatio::accumulate);
        zipper.addPositions(normalRatios, normalRatio::accumulate);
        zipper.run();

        return result;
    }

    @Override
    public void enter(GenomeRegion region) {
        baf.reset();
        tumorRatio.reset();
        normalRatio.reset();
    }

    @Override
    public void exit(GenomeRegion region) {
        double myTumorRatio = tumorRatio.meanRatio();
        double myNormalRatio = normalRatio.meanRatio();
        EnrichedCopyNumber copyNumber = ImmutableEnrichedCopyNumber.builder()
                .from(region)
                .mBAFCount(baf.count())
                .mBAF(baf.medianBaf())
                .tumorRatio(myTumorRatio)
                .normalRatio(myNormalRatio)
                .build();

        result.add(copyNumber);
    }

    private void variant(GermlineVariant variant) {
        final GermlineSampleData tumorData = variant.tumorData();
        if (tumorData == null || !GENO_TYPE.contains(variant.refData().genoType()) || variant.type() != VariantType.SNP
                || variant.refData().alleleFrequency() <= minRefAlleleFrequency
                || variant.refData().alleleFrequency() >= maxRefAlleleFrequency
                || variant.refData().combinedDepth() <= minCombinedDepth
                || variant.refData().combinedDepth() >= maxCombinedDepth) {
            return;
        }

        double standardBAH = tumorData.alleleFrequency();
        double modifiedBAF = 0.5 + Math.abs(standardBAH - 0.5);
        baf.accumulate(modifiedBAF);
    }

    private class BAFAccumulator {
        private int count;
        final private List<Double> bafs = Lists.newArrayList();

        private void reset() {
            count = 0;
            bafs.clear();
        }

        private void accumulate(double baf) {
            count++;
            bafs.add(baf);
        }

        private int count() {
            return count;
        }

        private double medianBaf() {
            if (count > 0) {
                Collections.sort(bafs);
                return bafs.size() % 2 == 0 ?
                        (bafs.get(count / 2) + bafs.get(count / 2 - 1)) / 2 :
                        bafs.get(count / 2);
            }
            return 0;
        }
    }

    private class RatioAccumulator {
        private double sumRatio;
        private int count;

        private void accumulate(FreecRatio ratio) {
            if (ratio.ratio() > -1) {
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

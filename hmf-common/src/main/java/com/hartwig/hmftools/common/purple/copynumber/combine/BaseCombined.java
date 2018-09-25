package com.hartwig.hmftools.common.purple.copynumber.combine;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ModifiableFittedRegion;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class BaseCombined implements GenomeRegion {

    private final ModifiableFittedRegion combined;
    private final List<FittedRegion> regions = Lists.newArrayList();

    public BaseCombined(final FittedRegion region) {
        this.combined = ModifiableFittedRegion.create().from(region);
        this.regions.add(region);
    }

    @NotNull
    @Override
    public String chromosome() {
        return combined.chromosome();
    }

    @Override
    public long start() {
        return combined.start();
    }

    @Override
    public long end() {
        return combined.end();
    }

    public FittedRegion region() {
        return combined;
    }

    public List<FittedRegion> regions() {
        return regions;
    }

    public void extend(final FittedRegion region) {
        combined.setStart(Math.min(combined.start(), region.start()));
        combined.setEnd(Math.max(combined.end(), region.end()));

        if (region.start() <= combined.start()) {
            regions.add(0, region);
            combined.setSupport(region.support());
            combined.setRatioSupport(region.ratioSupport());
        } else {
            regions.add(region);
        }
    }

    void applyBafWeights(final FittedRegion region, long currentWeight, long newWeight) {
        if (!Doubles.isZero(region.observedBAF())) {
            combined.setObservedBAF(weightedAverage(currentWeight, combined.observedBAF(), newWeight, region.observedBAF()));
        }

        if (!Doubles.isZero(region.tumorBAF())) {
            combined.setTumorBAF(weightedAverage(currentWeight, combined.tumorBAF(), newWeight, region.tumorBAF()));
        }

        combined.setBafCount(combined.bafCount() + region.bafCount());
    }

    void applyCopyNumberWeights(final FittedRegion region, long currentWeight, long newWeight) {

        if (!Doubles.isZero(region.tumorCopyNumber())) {
            combined.setTumorCopyNumber(weightedAverage(currentWeight, combined.tumorCopyNumber(), newWeight, region.tumorCopyNumber()));
        }

        if (!Doubles.isZero(region.refNormalisedCopyNumber())) {
            combined.setRefNormalisedCopyNumber(weightedAverage(currentWeight,
                    combined.refNormalisedCopyNumber(),
                    newWeight,
                    region.refNormalisedCopyNumber()));
        }

        combined.setDepthWindowCount(combined.depthWindowCount() + region.depthWindowCount());
    }

    double weightedAverage(long currentWeight, double currentValue, long newWeight, double newValue) {
        if (Doubles.isZero(currentValue)) {
            return newValue;
        }

        long totalWeight = currentWeight + newWeight;
        return (currentWeight * currentValue + newWeight * newValue) / totalWeight;
    }
}

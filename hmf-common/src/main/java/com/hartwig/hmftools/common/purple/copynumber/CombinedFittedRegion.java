package com.hartwig.hmftools.common.purple.copynumber;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ModifiableFittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegionStatus;

class CombinedFittedRegion {

    private final boolean bafWeighted;
    private ModifiableFittedRegion combined;
    private boolean modified = false;

    @Deprecated
    CombinedFittedRegion(final boolean bafWeighted, final FittedRegion region) {
        this(bafWeighted, region, region.status() != ObservedRegionStatus.SOMATIC);
    }

    CombinedFittedRegion(final boolean bafWeighted, final FittedRegion region, final boolean clearValues) {
        this.bafWeighted = bafWeighted;
        this.combined = ModifiableFittedRegion.create().from(region);
        if (clearValues) {
            clearValues();
        }
    }

    public void clearValues() {
        modified = true;
        combined.setRefNormalisedCopyNumber(0);
        combined.setObservedBAF(0);
        combined.setObservedTumorRatioCount(0);
        combined.setTumorCopyNumber(0);
        combined.setBafCount(0);
    }

    public boolean isModified() {
        return modified;
    }

    public void setModified() {
        this.modified = true;
    }

    FittedRegion region() {
        return combined;
    }

    @Deprecated
    void combine(final FittedRegion region) {
        combine(region, region.status() == ObservedRegionStatus.SOMATIC);
    }

    void combine(final FittedRegion region, boolean includeFittedValues) {
        modified = true;
        long currentBases = combined.bases();

        combined.setStart(Math.min(combined.start(), region.start()));
        combined.setEnd(Math.max(combined.end(), region.end()));

        if (region.start() <= combined.start()) {
            combined.setSupport(region.support());
            combined.setRatioSupport(region.ratioSupport());
        }

        if (includeFittedValues) {
            combined.setStatus(ObservedRegionStatus.SOMATIC); //TODO: Consider fixing this
            combined.setObservedTumorRatioCount(combined.observedTumorRatioCount() + region.observedTumorRatioCount());

            final long currentWeight;
            final long newWeight;
            if (bafWeighted && (combined.bafCount() > 0 || region.bafCount() > 0)) {
                currentWeight = combined.bafCount();
                newWeight = region.bafCount();
            } else {
                currentWeight = currentBases;
                newWeight = region.bases();
            }

            if (!Doubles.isZero(region.observedBAF())) {
                combined.setObservedBAF(newValue(currentWeight, combined.observedBAF(), newWeight, region.observedBAF()));
            }

            if (!Doubles.isZero(region.tumorCopyNumber())) {
                combined.setTumorCopyNumber(newValue(currentWeight, combined.tumorCopyNumber(), newWeight, region.tumorCopyNumber()));
            }

            if (!Doubles.isZero(region.refNormalisedCopyNumber())) {
                combined.setRefNormalisedCopyNumber(newValue(currentWeight,
                        combined.refNormalisedCopyNumber(),
                        newWeight,
                        region.refNormalisedCopyNumber()));
            }

            combined.setBafCount(combined.bafCount() + region.bafCount());
        }
    }

    private double newValue(long currentWeight, double currentValue, long newWeight, double newValue) {
        if (Doubles.isZero(currentValue)) {
            return newValue;
        }

        long totalWeight = currentWeight + newWeight;
        return (currentWeight * currentValue + newWeight * newValue) / totalWeight;
    }
}

package com.hartwig.hmftools.common.purple.copynumber;

import java.util.Optional;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantPloidy;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ModifiableFittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegionStatus;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

@SuppressWarnings("OptionalUsedAsFieldOrParameterType")
class CombinedRegion implements GenomeRegion {

    private final boolean bafWeighted;
    private ModifiableFittedRegion combined;
    private CombinedRegionMethod method = CombinedRegionMethod.NONE;
    private int unweightedCount = 0;

    @Deprecated
    CombinedRegion(final boolean bafWeighted, final FittedRegion region) {
        this(bafWeighted, region, region.status() != ObservedRegionStatus.SOMATIC);
    }

    CombinedRegion(final boolean bafWeighted, final FittedRegion region, final boolean clearValues) {
        this.bafWeighted = bafWeighted;
        this.combined = ModifiableFittedRegion.create().from(region);
        if (clearValues) {
            clearValues();
        }
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

    public double tumorCopyNumber() {
        return combined.tumorCopyNumber();
    }

    @Deprecated
    private void clearValues() {
        combined.setRefNormalisedCopyNumber(0);
        combined.setObservedBAF(0);
        combined.setObservedTumorRatioCount(0);
        combined.setTumorCopyNumber(0);
        combined.setBafCount(0);
    }

    public CombinedRegionMethod method() {
        return method;
    }

    public boolean isProcessed() {
        return method != CombinedRegionMethod.NONE;
    }

    public void setMethod(CombinedRegionMethod method) {
        this.method = method;
    }

    FittedRegion region() {
        return combined;
    }

    @Deprecated
    void combine(final FittedRegion region) {
        if (region.status() == ObservedRegionStatus.SOMATIC) {
            extendWithBAFWeightedAverage(region);
        } else {
            extend(region);
        }
    }

    void extend(final FittedRegion region) {
        combined.setStart(Math.min(combined.start(), region.start()));
        combined.setEnd(Math.max(combined.end(), region.end()));

        if (region.start() <= combined.start()) {
            combined.setSupport(region.support());
            combined.setRatioSupport(region.ratioSupport());
        }
    }

    void extendWithUnweightedAverage(final FittedRegion region) {
        extend(region);
        applyWeightedAverage(region, unweightedCount, 1);
        unweightedCount++;
    }

    void extendWithBAFWeightedAverage(final FittedRegion region) {
        long currentBases = combined.bases();

        extend(region);

        combined.setStatus(ObservedRegionStatus.SOMATIC);
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

        applyWeightedAverage(region, currentWeight, newWeight);
    }

    private void applyWeightedAverage(final FittedRegion region, long currentWeight, long newWeight) {

        if (!Doubles.isZero(region.observedBAF())) {
            combined.setObservedBAF(weightedAverage(currentWeight, combined.observedBAF(), newWeight, region.observedBAF()));
        }

        if (!Doubles.isZero(region.tumorBAF())) {
            combined.setTumorBAF(weightedAverage(currentWeight, combined.tumorBAF(), newWeight, region.tumorBAF()));
        }

        if (!Doubles.isZero(region.tumorCopyNumber())) {
            combined.setTumorCopyNumber(weightedAverage(currentWeight, combined.tumorCopyNumber(), newWeight, region.tumorCopyNumber()));
        }

        if (!Doubles.isZero(region.refNormalisedCopyNumber())) {
            combined.setRefNormalisedCopyNumber(weightedAverage(currentWeight,
                    combined.refNormalisedCopyNumber(),
                    newWeight,
                    region.refNormalisedCopyNumber()));
        }

        combined.setBafCount(combined.bafCount() + region.bafCount());
    }

    void inferCopyNumberFromStructuralVariants(final Optional<StructuralVariantPloidy> start, final Optional<StructuralVariantPloidy> end) {
        if (start.isPresent() || end.isPresent()) {
            setMethod(CombinedRegionMethod.STRUCTURAL_VARIANT);

            final double startWeight = start.map(StructuralVariantPloidy::impliedRightCopyNumberWeight).orElse(0d);
            final double startCopyNumber = start.map(StructuralVariantPloidy::impliedRightCopyNumber).orElse(0d);

            final double endWeight = end.map(StructuralVariantPloidy::impliedLeftCopyNumberWeight).orElse(0d);
            final double endCopyNumber = end.map(StructuralVariantPloidy::impliedLeftCopyNumber).orElse(0d);

            final double newCopyNumber = (startCopyNumber * startWeight + endCopyNumber * endWeight) / (startWeight + endWeight);
            combined.setTumorCopyNumber(newCopyNumber);
        }
    }

    private double weightedAverage(long currentWeight, double currentValue, long newWeight, double newValue) {
        if (Doubles.isZero(currentValue)) {
            return newValue;
        }

        long totalWeight = currentWeight + newWeight;
        return (currentWeight * currentValue + newWeight * newValue) / totalWeight;
    }
}

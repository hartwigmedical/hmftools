package com.hartwig.hmftools.common.purple.copynumber;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

class LowConfidenceCopyNumberBuilder extends BaseCopyNumberBuilder {

    private int bafCount;
    private double sumWeightedBAF;
    private double sumWeightedCopyNumber;
    private double sumWeightedRefNormalisedCopyNumber;

    LowConfidenceCopyNumberBuilder(double purity, @NotNull final FittedRegion fittedRegion) {
        super(purity, fittedRegion);
    }

    public int bafCount() {
        return bafCount;
    }

    public double averageObservedBAF() {
        return sumWeightedBAF / bases();
    }

    public double averageTumorCopyNumber() {
        return sumWeightedCopyNumber / bases();
    }

    @Override
    void extendRegion(@NotNull final FittedRegion value) {
        super.extendRegion(value);

        bafCount += value.bafCount();

        if (!Doubles.isZero(value.tumorCopyNumber())) {
            sumWeightedCopyNumber += value.tumorCopyNumber() * value.bases();
        }

        if (!Doubles.isZero(value.observedBAF())) {
            sumWeightedBAF += value.observedBAF() * value.bases();
        }
    }
}

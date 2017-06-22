package com.hartwig.hmftools.common.purple.copynumber;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

class LowConfidenceCopyNumberBuilder extends BaseCopyNumberBuilder {

    private int bafCount;
    private double sumWeightedBAF;
    private double sumWeightedCopyNumber;
    private double sumWeightedRefNormalisedCopyNumber;

    private long copyNumberBases;
    private long refNormalisedCopyNumberBases;
    private long bafBases;

    LowConfidenceCopyNumberBuilder(double purity, @NotNull final FittedRegion fittedRegion) {
        super(purity, fittedRegion);
    }

    @Override
    public int bafCount() {
        return bafCount;
    }

    @Override
    public double averageObservedBAF() {
        return bafBases == 0 ? 0 : sumWeightedBAF / bafBases;
    }

    @Override
    public double averageTumorCopyNumber() {
        return sumWeightedCopyNumber / copyNumberBases;
    }

    @Override
    public double averageRefNormalisedCopyNumber() {
        return sumWeightedRefNormalisedCopyNumber / refNormalisedCopyNumberBases;
    }

    @Override
    void extendRegion(@NotNull final FittedRegion value) {
        super.extendRegion(value);

        bafCount += value.bafCount();

        if (!Doubles.isZero(value.tumorCopyNumber())) {
            copyNumberBases += value.bases();
            sumWeightedCopyNumber += value.tumorCopyNumber() * value.bases();
        }

        if (!Doubles.isZero(value.observedBAF())) {
            refNormalisedCopyNumberBases += value.bases();
            sumWeightedBAF += value.observedBAF() * value.bases();
        }

        if (!Doubles.isZero(value.refNormalisedCopyNumber())) {
            bafBases += value.bases();
            sumWeightedRefNormalisedCopyNumber += value.refNormalisedCopyNumber() * value.bases();
        }
    }
}

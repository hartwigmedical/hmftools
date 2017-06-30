package com.hartwig.hmftools.common.purple.copynumber;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purity.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

class HighConfidenceCopyNumberBuilder extends BaseCopyNumberBuilder {

    private boolean weighWithBaf;

    private int totalBAFWeight;
    private int totalCopyNumberWeight;

    private double sumWeightedBAF;
    private double sumWeightedCopyNumber;
    private double sumWeightedRefNormalisedCopyNumber;

    HighConfidenceCopyNumberBuilder(@NotNull final PurityAdjuster purityAdjuster, @NotNull final FittedRegion fittedRegion) {
        super(purityAdjuster, fittedRegion);
    }

    @Override
    public int bafCount() {
        return weighWithBaf ? totalBAFWeight : 0;
    }

    @Override
    public double averageObservedBAF() {
        return totalBAFWeight == 0 ? 0 : sumWeightedBAF / totalBAFWeight;
    }

    @Override
    public double averageTumorCopyNumber() {
        return totalCopyNumberWeight == 0 ? 0 : sumWeightedCopyNumber / totalCopyNumberWeight;
    }

    @Override
    public double averageRefNormalisedCopyNumber() {
        return totalCopyNumberWeight == 0 ? 0 : sumWeightedRefNormalisedCopyNumber / totalCopyNumberWeight;
    }

    public void extendRegion(@NotNull final FittedRegion value) {
        super.extendRegion(value);

        if (value.bafCount() > 0) {

            if (!weighWithBaf) {
                resetAverage();
                weighWithBaf = true;
            }

            long weight = value.bafCount();
            averageInBaf(weight, value.observedBAF());
            averageInCopyNumber(weight, value.tumorCopyNumber(), value.refNormalisedCopyNumber());

        } else if (!weighWithBaf) {
            long weight = Math.max(1, value.bases() / 1000);
            averageInCopyNumber(weight, value.tumorCopyNumber(), value.refNormalisedCopyNumber());
        }
    }

    private void averageInBaf(long weight, double baf) {
        totalBAFWeight += weight;
        sumWeightedBAF += baf * weight;
    }

    private void averageInCopyNumber(long weight, double copyNumber, double refNormalisedCopyNumber) {
        if (!Doubles.isZero(copyNumber)) {
            totalCopyNumberWeight += weight;
            sumWeightedCopyNumber += copyNumber * weight;
            sumWeightedRefNormalisedCopyNumber += refNormalisedCopyNumber * weight;
        }
    }

    private void resetAverage() {
        totalBAFWeight = 0;
        totalCopyNumberWeight = 0;

        sumWeightedBAF = 0;
        sumWeightedCopyNumber = 0;
        sumWeightedRefNormalisedCopyNumber = 0;
    }
}

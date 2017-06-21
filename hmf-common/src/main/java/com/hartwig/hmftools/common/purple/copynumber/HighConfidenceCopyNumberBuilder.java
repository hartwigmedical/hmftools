package com.hartwig.hmftools.common.purple.copynumber;

import static com.hartwig.hmftools.common.numeric.Doubles.lessOrEqual;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purity.PurityAdjustment;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.FittedRegionFactory;

import org.jetbrains.annotations.NotNull;

class HighConfidenceCopyNumberBuilder {

    private final double purity;
    private final String chromosome;
    private long start = 1;
    private long end;

    private boolean weighWithBaf;

    private int totalBAFWeight;
    private int totalCopyNumberWeight;

    private double sumWeightedBAF;
    private double sumWeightedCopyNumber;
    private double sumWeightedRefNormalisedCopyNumber;

    HighConfidenceCopyNumberBuilder(double purity, @NotNull final FittedRegion fittedRegion) {
        this.purity = purity;
        this.chromosome = fittedRegion.chromosome();
        this.start = fittedRegion.start();
        extendRegion(fittedRegion);
    }

    public String chromosome() {
        return chromosome;
    }

    private int bafCount() {
        return weighWithBaf ? totalBAFWeight : 0;
    }

    double averageObservedBAF() {
        return totalBAFWeight == 0 ? 0 : sumWeightedBAF / totalBAFWeight;
    }

    double averageTumorCopyNumber() {
        return totalCopyNumberWeight == 0 ? 0 : sumWeightedCopyNumber / totalCopyNumberWeight;
    }

    double averageRefNormalisedCopyNumber() {
        return totalCopyNumberWeight == 0 ? 0 : sumWeightedRefNormalisedCopyNumber / totalCopyNumberWeight;
    }

    void extendRegion(@NotNull final FittedRegion value) {
        assert (chromosome.equals(value.chromosome())) : "Regions cannot be extended between chromosomes";

        start = Math.min(value.start(), start);
        end = Math.max(value.end(), end);

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

    @NotNull
    public PurpleCopyNumber build() {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .bafCount(bafCount())
                .averageObservedBAF(averageObservedBAF())
                .averageActualBAF(purityAdjustedBAF(purity, averageTumorCopyNumber(), averageObservedBAF()))
                .averageTumorCopyNumber(averageTumorCopyNumber())
                .build();
    }

    @VisibleForTesting
    static double purityAdjustedBAF(final double purity, final double copyNumber, final double observedBAF) {
        if (Doubles.isZero(copyNumber)) {
            return observedBAF;
        }

        double adjustedObservedBAF =
                isEven(copyNumber) && lessOrEqual(observedBAF, FittedRegionFactory.NORMAL_BAF) ? 0.5 : observedBAF;
        return PurityAdjustment.purityAdjustedBAF(purity, copyNumber, adjustedObservedBAF);
    }

    @VisibleForTesting
    static boolean isEven(double copyNumber) {

        double decimal = copyNumber % 1d;
        double wholeNumber = copyNumber - decimal;

        return (wholeNumber % 2 == 0 && Doubles.lessOrEqual(decimal, 0.25)) || (wholeNumber % 2 != 0
                && Doubles.greaterOrEqual(decimal, 0.75));
    }
}

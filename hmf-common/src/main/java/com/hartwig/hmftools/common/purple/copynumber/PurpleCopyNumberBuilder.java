package com.hartwig.hmftools.common.purple.copynumber;

import static com.hartwig.hmftools.common.numeric.Doubles.lessOrEqual;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purity.PurityAdjustment;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.FittedRegionFactory;

import org.jetbrains.annotations.NotNull;

class PurpleCopyNumberBuilder {

    private final double purity;
    private final String chromosome;
    private long start = 1;
    private long end;

    private boolean weighWithBaf;
    private int totalWeight;
    private double sumWeightedBAF;
    private double sumWeightedCopyNumber;
    private double sumWeightedRefNormalisedCopyNumber;

    PurpleCopyNumberBuilder(double purity, @NotNull final FittedRegion fittedRegion) {
        this.purity = purity;
        this.chromosome = fittedRegion.chromosome();
        this.start = fittedRegion.start();
        extendRegion(fittedRegion);
    }

    public String chromosome() {
        return chromosome;
    }

    private int bafCount() {
        return weighWithBaf ? totalWeight : 0;
    }

    double averageObservedBAF() {
        return totalWeight == 0 ? 0 : sumWeightedBAF / totalWeight;
    }

    double averageTumorCopyNumber() {
        return totalWeight == 0 ? 0 : sumWeightedCopyNumber / totalWeight;
    }

    double averageRefNormalisedCopyNumber() {
        return totalWeight == 0 ? 0 : sumWeightedRefNormalisedCopyNumber / totalWeight;
    }

    void extendRegion(@NotNull final FittedRegion value) {
        assert (chromosome.equals(value.chromosome())) : "Regions cannot be extended between chromosomes";

        start = Math.min(value.start(), start);
        end = Math.max(value.end(), end);

        double ratio = value.tumorCopyNumber();
        double baf = value.observedBAF();

        if (value.bafCount() > 0) {

            if (!weighWithBaf) {
                resetAverage();
                weighWithBaf = true;
            }

            long weight = value.bafCount();
            totalWeight += weight;
            sumWeightedBAF += baf * weight;
            sumWeightedCopyNumber += ratio * weight;
            sumWeightedRefNormalisedCopyNumber += value.refNormalisedCopyNumber() * weight;

        } else if (!weighWithBaf && !Doubles.isZero(ratio)) {

            long weight = Math.max(1, value.bases() / 1000);
            totalWeight += weight;
            sumWeightedCopyNumber += ratio * weight;
            sumWeightedRefNormalisedCopyNumber += value.refNormalisedCopyNumber() * weight;
        }
    }

    private void resetAverage() {
        totalWeight = 0;
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
        double adjustedObservedBAF = isEven(copyNumber) && lessOrEqual(observedBAF, FittedRegionFactory.NORMAL_BAF)
                ? 0.5
                : observedBAF;
        return PurityAdjustment.purityAdjustedFrequency(purity, copyNumber, adjustedObservedBAF, 0.5);
    }

    @VisibleForTesting
    static boolean isEven(double copyNumber) {

        double decimal = copyNumber % 1d;
        double wholeNumber = copyNumber - decimal;

        return (wholeNumber % 2 == 0 && Doubles.lessOrEqual(decimal, 0.25)) || (wholeNumber % 2 != 0
                && Doubles.greaterOrEqual(decimal, 0.75));
    }
}

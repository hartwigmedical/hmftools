package com.hartwig.hmftools.common.purple.copynumber;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.FittedRegion;

import org.jetbrains.annotations.NotNull;

class PurpleCopyNumberBuilder {

    private final String chromosome;
    private long start = 1;
    private long end;

    private boolean weighWithBaf;
    private int totalWeight;
    private double sumWeightedBAF;
    private double sumWeightedPurityAdjustedBAF;
    private double sumWeightedCopyNumber;
    private double sumWeightedRefNormalisedCopyNumber;

    PurpleCopyNumberBuilder(@NotNull final FittedRegion fittedRegion) {
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

    private double averagePurityAdjustedBAF() {
        return totalWeight == 0 ? 0 : sumWeightedPurityAdjustedBAF / totalWeight;
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
            sumWeightedPurityAdjustedBAF += value.purityAdjustedBAF() * weight;
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
        sumWeightedPurityAdjustedBAF = 0;
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
                .averageActualBAF(averagePurityAdjustedBAF())
                .averageTumorCopyNumber(averageTumorCopyNumber())
                .build();
    }
}

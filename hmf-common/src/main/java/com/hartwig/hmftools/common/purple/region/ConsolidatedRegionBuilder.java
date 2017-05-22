package com.hartwig.hmftools.common.purple.region;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.FittedCopyNumber;

class ConsolidatedRegionBuilder {

    private final String chromosome;
    private long start = 1;
    private long end;

    private boolean weighWithBaf;
    private long totalWeight;
    private double sumWeightedBAF;
    private double sumWeightedRatioOfRatios;

    ConsolidatedRegionBuilder(FittedCopyNumber copyNumber) {
        this.chromosome = copyNumber.chromosome();
        this.start = copyNumber.start();
        extendRegion(copyNumber);
    }

    public String chromosome() {
        return chromosome;
    }

    double averageBAF() {
        return totalWeight == 0 ? 0 : sumWeightedBAF / totalWeight;
    }

    double averageRatioOfRatios() {
        return totalWeight == 0 ? 0 : sumWeightedRatioOfRatios / totalWeight;
    }

    void extendRegion(FittedCopyNumber value) {
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
            sumWeightedRatioOfRatios += ratio * weight;

        } else if (!weighWithBaf && !Doubles.isZero(ratio)) {

            long weight = Math.max(1, value.bases()/1000);
            totalWeight += weight;
            sumWeightedRatioOfRatios += ratio * weight;
        }
    }

    private void resetAverage() {
        totalWeight = 0;
        sumWeightedBAF = 0;
        sumWeightedRatioOfRatios = 0;
    }

    public ConsolidatedRegion build() {
        return ImmutableConsolidatedRegion.builder().chromosome(chromosome).start(start).end(end).averageBAF(
                averageBAF()).averageRatioOfRatios(averageRatioOfRatios()).build();
    }
}

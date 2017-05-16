package com.hartwig.hmftools.common.purple;

class ConsolidatedRegionBuilder {

    private final String chromosome;
    private long start = 1;
    private long end;

    private long weight;
    private double sumWeightedBAF;
    private double sumWeightedRatioOfRatios;

    ConsolidatedRegionBuilder(final String chromosome, long start) {
        this.chromosome = chromosome;
        this.start = start;
    }

    ConsolidatedRegionBuilder(FittedCopyNumber copyNumber) {
        this.chromosome = copyNumber.chromosome();
        this.start = copyNumber.start();
        extendRegion(copyNumber);
    }


    public double averageBAF() {
        return sumWeightedBAF / weight;
    }

    public double averageRatioOfRatios() {
        return sumWeightedRatioOfRatios / weight;
    }

    public void extendRegion(FittedCopyNumber value) {
        assert (chromosome.equals(value.chromosome())) : "Regions cannot be extended between chromosomes";

        start = Math.min(value.start(), start);
        end = Math.max(value.end(), end);

        if (value.bafCount() > 0) {
            weight += value.bafCount();
            sumWeightedBAF += value.actualBAF() * value.bafCount();
            sumWeightedRatioOfRatios += value.ratioOfRatios() * value.bafCount();
        }
    }

    public ConsolidatedRegion build() {
        return ImmutableConsolidatedRegion.builder()
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .averageBAF(averageBAF())
                .averageRatioOfRatios(averageRatioOfRatios())
                .build();
    }
}

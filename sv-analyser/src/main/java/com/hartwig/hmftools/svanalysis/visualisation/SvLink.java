package com.hartwig.hmftools.svanalysis.visualisation;

public class SvLink {


    private final String startChromosome;
    private final long startPosition;
    private final String endChromosome;
    private final long endPosition;

    public SvLink(final String startChromosome, final long startPosition, final String endChromosome, final long endPosition) {
        this.startChromosome = startChromosome;
        this.startPosition = startPosition;
        this.endChromosome = endChromosome;
        this.endPosition = endPosition;
    }

    public String startChromosome() {
        return startChromosome;
    }

    public long startPosition() {
        return startPosition;
    }

    public String endChromosome() {
        return endChromosome;
    }

    public long endPosition() {
        return endPosition;
    }

    @Override
    public String toString() {
        return "SvLink{" + "startChromosome='" + startChromosome + '\'' + ", startPosition=" + startPosition + ", endChromosome='"
                + endChromosome + '\'' + ", endPosition=" + endPosition + '}';
    }
}

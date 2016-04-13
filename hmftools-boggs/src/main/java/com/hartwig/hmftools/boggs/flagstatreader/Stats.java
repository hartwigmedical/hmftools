package com.hartwig.hmftools.boggs.flagstatreader;

public class Stats {
    private final int total;
    private final int secondary;
    private final int supplementary;
    private final int duplicates;
    private final int mapped;
    private final int pairedInSequencing;
    private final int read1;
    private final int read2;
    private final int properlyPaired;
    private final int itselfAndMateMapped;
    private final int singletons;
    private final int mateMappedToDifferentChr;
    private final int mateMappedToDifferentChrMapQ5;

    Stats(int total, int secondary, int supplementary, int duplicates, int mapped,
                 int pairedInSequencing, int read1, int read2, int properlyPaired, int itselfAndMateMapped,
                 int singletons, int mateMappedToDifferentChr, int mateMappedToDifferentChrMapQ5) {
        this.total = total;
        this.secondary = secondary;
        this.supplementary = supplementary;
        this.duplicates = duplicates;
        this.mapped = mapped;
        this.pairedInSequencing = pairedInSequencing;
        this.read1 = read1;
        this.read2 = read2;
        this.properlyPaired = properlyPaired;
        this.itselfAndMateMapped = itselfAndMateMapped;
        this.singletons = singletons;
        this.mateMappedToDifferentChr = mateMappedToDifferentChr;
        this.mateMappedToDifferentChrMapQ5 = mateMappedToDifferentChrMapQ5;
    }

    public int getTotal() {
        return total;
    }

    public int getSecondary() {
        return secondary;
    }

    public int getSupplementary() {
        return supplementary;
    }

    public int getDuplicates() {
        return duplicates;
    }

    public int getMapped() {
        return mapped;
    }

    public int getPairedInSequencing() {
        return pairedInSequencing;
    }

    public int getRead1() {
        return read1;
    }

    public int getRead2() {
        return read2;
    }

    public int getProperlyPaired() {
        return properlyPaired;
    }

    public int getItselfAndMateMapped() {
        return itselfAndMateMapped;
    }

    public int getSingletons() {
        return singletons;
    }

    public int getMateMappedToDifferentChr() {
        return mateMappedToDifferentChr;
    }

    public int getMateMappedToDifferentChrMapQ5() {
        return mateMappedToDifferentChrMapQ5;
    }
}

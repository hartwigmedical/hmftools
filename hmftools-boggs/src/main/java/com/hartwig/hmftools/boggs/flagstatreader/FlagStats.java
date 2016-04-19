package com.hartwig.hmftools.boggs.flagstatreader;

public class FlagStats {
    private final long total;
    private final long secondary;
    private final long supplementary;
    private final long duplicates;
    private final long mapped;
    private final long pairedInSequencing;
    private final long read1;
    private final long read2;
    private final long properlyPaired;
    private final long itselfAndMateMapped;
    private final long singletons;
    private final long mateMappedToDifferentChr;
    private final long mateMappedToDifferentChrMapQ5;

    FlagStats(long total, long secondary, long supplementary, long duplicates, long mapped, long pairedInSequencing,
                     long read1, long read2, long properlyPaired, long itselfAndMateMapped, long singletons,
                     long mateMappedToDifferentChr, long mateMappedToDifferentChrMapQ5) {
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

    public long total() {
        return total;
    }

    public long secondary() {
        return secondary;
    }

    public long supplementary() {
        return supplementary;
    }

    public long duplicates() {
        return duplicates;
    }

    public long mapped() {
        return mapped;
    }

    public long pairedInSequencing() {
        return pairedInSequencing;
    }

    public long read1() {
        return read1;
    }

    public long read2() {
        return read2;
    }

    public long properlyPaired() {
        return properlyPaired;
    }

    public long itselfAndMateMapped() {
        return itselfAndMateMapped;
    }

    public long singletons() {
        return singletons;
    }

    public long mateMappedToDifferentChr() {
        return mateMappedToDifferentChr;
    }

    public long mateMappedToDifferentChrMapQ5() {
        return mateMappedToDifferentChrMapQ5;
    }

    @Override
    public String toString() {
        return "FlagStats{" +
                "total=" + total +
                ", secondary=" + secondary +
                ", supplementary=" + supplementary +
                ", duplicates=" + duplicates +
                ", mapped=" + mapped +
                ", pairedInSequencing=" + pairedInSequencing +
                ", read1=" + read1 +
                ", read2=" + read2 +
                ", properlyPaired=" + properlyPaired +
                ", itselfAndMateMapped=" + itselfAndMateMapped +
                ", singletons=" + singletons +
                ", mateMappedToDifferentChr=" + mateMappedToDifferentChr +
                ", mateMappedToDifferentChrMapQ5=" + mateMappedToDifferentChrMapQ5 +
                '}';
    }
}

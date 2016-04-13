package com.hartwig.hmftools.boggs.flagstatreader;

import org.jetbrains.annotations.NotNull;

public class StatsBuilder {
    private long total;
    private long secondary;
    private long supplementary;
    private long duplicates;
    private long mapped;
    private long pairedInSequencing;
    private long read1;
    private long read2;
    private long properlyPaired;
    private long itselfAndMateMapped;
    private long singletons;
    private long mateMappedToDifferentChr;
    private long mateMappedToDifferentChrMapQ5;

    @NotNull
    public StatsBuilder setTotal(long total) {
        this.total = total;
        return this;
    }

    @NotNull
    public StatsBuilder setSecondary(long secondary) {
        this.secondary = secondary;
        return this;
    }

    @NotNull
    public StatsBuilder setSupplementary(long supplementary) {
        this.supplementary = supplementary;
        return this;
    }

    @NotNull
    public StatsBuilder setDuplicates(long duplicates) {
        this.duplicates = duplicates;
        return this;
    }

    @NotNull
    public StatsBuilder setMapped(long mapped) {
        this.mapped = mapped;
        return this;
    }

    @NotNull
    public StatsBuilder setPairedInSequencing(long pairedInSequencing) {
        this.pairedInSequencing = pairedInSequencing;
        return this;
    }

    @NotNull
    public StatsBuilder setRead1(long read1) {
        this.read1 = read1;
        return this;
    }

    @NotNull
    public StatsBuilder setRead2(long read2) {
        this.read2 = read2;
        return this;
    }

    @NotNull
    public StatsBuilder setProperlyPaired(long properlyPaired) {
        this.properlyPaired = properlyPaired;
        return this;
    }

    @NotNull
    public StatsBuilder setItselfAndMateMapped(long itselfAndMateMapped) {
        this.itselfAndMateMapped = itselfAndMateMapped;
        return this;
    }

    @NotNull
    public StatsBuilder setSingletons(long singletons) {
        this.singletons = singletons;
        return this;
    }

    @NotNull
    public StatsBuilder setMateMappedToDifferentChr(long mateMappedToDifferentChr) {
        this.mateMappedToDifferentChr = mateMappedToDifferentChr;
        return this;
    }

    @NotNull
    public StatsBuilder setMateMappedToDifferentChrMapQ5(long mateMappedToDifferentChrMapQ5) {
        this.mateMappedToDifferentChrMapQ5 = mateMappedToDifferentChrMapQ5;
        return this;
    }

    @NotNull
    public Stats build() {
        return new Stats(total, secondary, supplementary, duplicates, mapped, pairedInSequencing ,read1, read2,
                properlyPaired, itselfAndMateMapped, singletons, mateMappedToDifferentChr, mateMappedToDifferentChrMapQ5);
    }
}

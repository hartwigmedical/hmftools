package com.hartwig.hmftools.boggs.flagstatreader;

import org.jetbrains.annotations.NotNull;

class StatsBuilder {
    private int total;
    private int secondary;
    private int supplementary;
    private int duplicates;
    private int mapped;
    private int pairedInSequencing;
    private int read1;
    private int read2;
    private int properlyPaired;
    private int itselfAndMateMapped;
    private int singletons;
    private int mateMappedToDifferentChr;
    private int mateMappedToDifferentChrMapQ5;

    @NotNull
    public StatsBuilder setTotal(int total) {
        this.total = total;
        return this;
    }

    @NotNull
    public StatsBuilder setSecondary(int secondary) {
        this.secondary = secondary;
        return this;
    }

    @NotNull
    public StatsBuilder setSupplementary(int supplementary) {
        this.supplementary = supplementary;
        return this;
    }

    @NotNull
    public StatsBuilder setDuplicates(int duplicates) {
        this.duplicates = duplicates;
        return this;
    }

    @NotNull
    public StatsBuilder setMapped(int mapped) {
        this.mapped = mapped;
        return this;
    }

    @NotNull
    public StatsBuilder setPairedInSequencing(int pairedInSequencing) {
        this.pairedInSequencing = pairedInSequencing;
        return this;
    }

    @NotNull
    public StatsBuilder setRead1(int read1) {
        this.read1 = read1;
        return this;
    }

    @NotNull
    public StatsBuilder setRead2(int read2) {
        this.read2 = read2;
        return this;
    }

    @NotNull
    public StatsBuilder setProperlyPaired(int properlyPaired) {
        this.properlyPaired = properlyPaired;
        return this;
    }

    @NotNull
    public StatsBuilder setItselfAndMateMapped(int itselfAndMateMapped) {
        this.itselfAndMateMapped = itselfAndMateMapped;
        return this;
    }

    @NotNull
    public StatsBuilder setSingletons(int singletons) {
        this.singletons = singletons;
        return this;
    }

    @NotNull
    public StatsBuilder setMateMappedToDifferentChr(int mateMappedToDifferentChr) {
        this.mateMappedToDifferentChr = mateMappedToDifferentChr;
        return this;
    }

    @NotNull
    public StatsBuilder setMateMappedToDifferentChrMapQ5(int mateMappedToDifferentChrMapQ5) {
        this.mateMappedToDifferentChrMapQ5 = mateMappedToDifferentChrMapQ5;
        return this;
    }

    @NotNull
    public Stats build() {
        return new Stats(total, secondary, supplementary, duplicates, mapped, pairedInSequencing ,read1, read2,
                properlyPaired, itselfAndMateMapped, singletons, mateMappedToDifferentChr, mateMappedToDifferentChrMapQ5);
    }
}

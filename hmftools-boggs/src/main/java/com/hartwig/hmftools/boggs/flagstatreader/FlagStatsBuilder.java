package com.hartwig.hmftools.boggs.flagstatreader;

import org.jetbrains.annotations.NotNull;

class FlagStatsBuilder {
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
    public FlagStatsBuilder setTotal(long total) {
        this.total = total;
        return this;
    }

    @NotNull
    public FlagStatsBuilder setSecondary(long secondary) {
        this.secondary = secondary;
        return this;
    }

    @NotNull
    public FlagStatsBuilder setSupplementary(long supplementary) {
        this.supplementary = supplementary;
        return this;
    }

    @NotNull
    public FlagStatsBuilder setDuplicates(long duplicates) {
        this.duplicates = duplicates;
        return this;
    }

    @NotNull
    public FlagStatsBuilder setMapped(long mapped) {
        this.mapped = mapped;
        return this;
    }

    @NotNull
    public FlagStatsBuilder setPairedInSequencing(long pairedInSequencing) {
        this.pairedInSequencing = pairedInSequencing;
        return this;
    }

    @NotNull
    public FlagStatsBuilder setRead1(long read1) {
        this.read1 = read1;
        return this;
    }

    @NotNull
    public FlagStatsBuilder setRead2(long read2) {
        this.read2 = read2;
        return this;
    }

    @NotNull
    public FlagStatsBuilder setProperlyPaired(long properlyPaired) {
        this.properlyPaired = properlyPaired;
        return this;
    }

    @NotNull
    public FlagStatsBuilder setItselfAndMateMapped(long itselfAndMateMapped) {
        this.itselfAndMateMapped = itselfAndMateMapped;
        return this;
    }

    @NotNull
    public FlagStatsBuilder setSingletons(long singletons) {
        this.singletons = singletons;
        return this;
    }

    @NotNull
    public FlagStatsBuilder setMateMappedToDifferentChr(long mateMappedToDifferentChr) {
        this.mateMappedToDifferentChr = mateMappedToDifferentChr;
        return this;
    }

    @NotNull
    public FlagStatsBuilder setMateMappedToDifferentChrMapQ5(long mateMappedToDifferentChrMapQ5) {
        this.mateMappedToDifferentChrMapQ5 = mateMappedToDifferentChrMapQ5;
        return this;
    }

    @NotNull
    public FlagStats build() {
        return new FlagStats(total, secondary, supplementary, duplicates, mapped, pairedInSequencing ,read1, read2,
                properlyPaired, itselfAndMateMapped, singletons, mateMappedToDifferentChr, mateMappedToDifferentChrMapQ5);
    }
}

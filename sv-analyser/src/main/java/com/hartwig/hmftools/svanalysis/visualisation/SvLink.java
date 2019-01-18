package com.hartwig.hmftools.svanalysis.visualisation;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

class SvLink implements GenomeRegion {

    private final GenomeRegion region;
    private final int value;

    SvLink(final GenomeRegion region, final int value) {
        this.region = region;
        this.value = value;
    }

    @NotNull
    @Override
    public String chromosome() {
        return region.chromosome();
    }

    @Override
    public long start() {
        return region.start();
    }

    @Override
    public long end() {
        return region.end();
    }

    public int value() {
        return value;
    }

}

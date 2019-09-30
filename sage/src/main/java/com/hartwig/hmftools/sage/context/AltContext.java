package com.hartwig.hmftools.sage.context;

import com.hartwig.hmftools.common.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

public class AltContext implements VariantHotspot {

    private final RefContext refContext;
    private final String alt;

    private int altQuality;
    private int altReads;

    public AltContext(final RefContext refContext, final String alt) {
        this.refContext = refContext;
        this.alt = alt;
    }

    @NotNull
    @Override
    public String ref() {
        return refContext.ref();
    }

    @NotNull
    @Override
    public String alt() {
        return alt;
    }

    @NotNull
    @Override
    public String chromosome() {
        return refContext.chromosome();
    }

    @Override
    public long position() {
        return refContext.position();
    }

    public int altQuality() {
        return altQuality;
    }

    public int altReads() {
        return altReads;
    }

    public int refQuality() {
        return refContext.refQuality();
    }

    public int refReads() {
        return refContext.refReads();
    }

    public int readDepth() {
        return refContext.readDepth();
    }

    public int subprimeReadDepth() {
        return refContext.subprimeReadDepth();
    }

}

package com.hartwig.hmftools.sage.candidate;

import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.read.ReadContext;

import org.jetbrains.annotations.NotNull;

public class Candidate {

    private final VariantHotspot variant;

    private int maxDepth;
    private int readContextSupport;
    private ReadContext readContext;

    public Candidate(final VariantHotspot variant) {
        this.variant = variant;
    }

    public Candidate(final AltContext altContext) {
        this.variant = ImmutableVariantHotspotImpl.builder().from(altContext).build();
        this.maxDepth = altContext.rawDepth();
        this.readContext = altContext.primaryReadContext().readContext();
        this.readContextSupport = altContext.primaryReadContext().altSupport();
    }

    public void update(final AltContext altContext) {
        int altContextSupport = altContext.primaryReadContext().altSupport();
        if (altContextSupport > readContextSupport) {
            readContextSupport = altContextSupport;
            readContext = altContext.primaryReadContext().readContext();
        }
        maxDepth = Math.max(maxDepth, altContext.rawDepth());
    }

    @NotNull
    public VariantHotspot variant() {
        return variant;
    }

    public int maxDepth() {
        return maxDepth;
    }

    @NotNull
    public ReadContext readContext() {
        return readContext;
    }
}

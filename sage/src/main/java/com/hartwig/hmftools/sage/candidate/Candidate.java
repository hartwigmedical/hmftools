package com.hartwig.hmftools.sage.candidate;

import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.read.ReadContext;

import org.jetbrains.annotations.NotNull;

public class Candidate {

    private final VariantHotspot hotspot;

    private int maxDepth;
    private int readContextSupport;
    private ReadContext readContext;

    public Candidate(final VariantHotspot hotspot) {
        this.hotspot = hotspot;
    }

    public Candidate(final AltContext altContext) {
        this.hotspot = ImmutableVariantHotspotImpl.builder().from(altContext).build();
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
    public VariantHotspot hotspot() {
        return hotspot;
    }

    public int maxDepth() {
        return maxDepth;
    }

    @NotNull
    public ReadContext readContext() {
        return readContext;
    }
}

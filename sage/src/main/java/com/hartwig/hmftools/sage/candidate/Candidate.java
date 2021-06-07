package com.hartwig.hmftools.sage.candidate;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.variant.SageVariantTier;

import org.jetbrains.annotations.NotNull;

public class Candidate implements GenomePosition
{

    private final SageVariantTier tier;
    private final VariantHotspot variant;

    private int maxDepth;
    private int minNumberOfEvents;
    private int readContextSupport;
    private ReadContext readContext;

    public Candidate(final SageVariantTier tier, final VariantHotspot variant, final ReadContext readContext, int maxDepth,
            int minNumberOfEvents)
    {
        this.tier = tier;
        this.variant = variant;
        this.readContext = readContext;
        this.maxDepth = maxDepth;
        this.minNumberOfEvents = minNumberOfEvents;
    }

    public Candidate(final SageVariantTier tier, final AltContext altContext)
    {
        this.tier = tier;
        this.variant = ImmutableVariantHotspotImpl.builder().from(altContext).build();
        this.maxDepth = altContext.rawDepth();
        this.readContext = altContext.readContext();
        this.readContextSupport = altContext.readContextSupport();
        this.minNumberOfEvents = altContext.minNumberOfEvents();
    }

    public void update(final AltContext altContext)
    {
        int altContextSupport = altContext.readContextSupport();
        if(altContextSupport > readContextSupport)
        {
            readContextSupport = altContextSupport;
            readContext = altContext.readContext();
            minNumberOfEvents = Math.min(minNumberOfEvents, altContext.minNumberOfEvents());
        }
        maxDepth = Math.max(maxDepth, altContext.rawDepth());
    }

    @NotNull
    public SageVariantTier tier()
    {
        return tier;
    }

    @NotNull
    public VariantHotspot variant()
    {
        return variant;
    }

    public int maxReadDepth()
    {
        return maxDepth;
    }

    @NotNull
    public ReadContext readContext()
    {
        return readContext;
    }

    public int minNumberOfEvents()
    {
        return minNumberOfEvents;
    }

    @NotNull
    @Override
    public String chromosome()
    {
        return variant.chromosome();
    }

    @Override
    public long position()
    {
        return variant.position();
    }
}

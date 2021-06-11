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
    private final SageVariantTier mTier;
    private final VariantHotspot mVariant;

    private int mMaxDepth;
    private int mMinNumberOfEvents;
    private int mReadContextSupport;
    private ReadContext mReadContext;

    public Candidate(final SageVariantTier tier, final VariantHotspot variant, final ReadContext readContext, int maxDepth,
            int minNumberOfEvents)
    {
        mTier = tier;
        mVariant = variant;
        mReadContext = readContext;
        mMaxDepth = maxDepth;
        mMinNumberOfEvents = minNumberOfEvents;
    }

    public Candidate(final SageVariantTier tier, final AltContext altContext)
    {
        mTier = tier;
        mVariant = ImmutableVariantHotspotImpl.builder().from(altContext).build();
        mMaxDepth = altContext.rawDepth();
        mReadContext = altContext.readContext();
        mReadContextSupport = altContext.readContextSupport();
        mMinNumberOfEvents = altContext.minNumberOfEvents();
    }

    public void update(final AltContext altContext)
    {
        int altContextSupport = altContext.readContextSupport();
        if(altContextSupport > mReadContextSupport)
        {
            mReadContextSupport = altContextSupport;
            mReadContext = altContext.readContext();
            mMinNumberOfEvents = Math.min(mMinNumberOfEvents, altContext.minNumberOfEvents());
        }
        mMaxDepth = Math.max(mMaxDepth, altContext.rawDepth());
    }

    @NotNull
    public SageVariantTier tier()
    {
        return mTier;
    }

    @NotNull
    public VariantHotspot variant()
    {
        return mVariant;
    }

    public int maxReadDepth()
    {
        return mMaxDepth;
    }

    @NotNull
    public ReadContext readContext()
    {
        return mReadContext;
    }

    public int minNumberOfEvents()
    {
        return mMinNumberOfEvents;
    }

    @NotNull
    @Override
    public String chromosome()
    {
        return mVariant.chromosome();
    }

    @Override
    public long position()
    {
        return mVariant.position();
    }
}

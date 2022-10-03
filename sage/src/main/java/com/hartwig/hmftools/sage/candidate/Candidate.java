package com.hartwig.hmftools.sage.candidate;

import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.common.VariantTier;

public class Candidate
{
    private final VariantTier mTier;
    private final VariantHotspot mVariant;

    private int mMinNumberOfEvents;
    private int mReadContextSupport;
    private ReadContext mReadContext;

    public Candidate(
            final VariantTier tier, final VariantHotspot variant, final ReadContext readContext, int minNumberOfEvents, int readContextSupport)
    {
        mTier = tier;
        mVariant = variant;
        mReadContext = readContext;
        mMinNumberOfEvents = minNumberOfEvents;
        mReadContextSupport = readContextSupport;
    }

    public static Candidate fromAltContext(final VariantTier tier, final AltContext altContext)
    {
        return new Candidate(
                tier, ImmutableVariantHotspotImpl.builder().from(altContext).build(), altContext.readContext(),
                altContext.minNumberOfEvents(), altContext.readContextSupport());
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
    }

    public VariantTier tier() { return mTier; }
    public VariantHotspot variant() { return mVariant; }

    public ReadContext readContext() { return mReadContext; }

    public int minNumberOfEvents() { return mMinNumberOfEvents; }

    public String chromosome() { return mVariant.chromosome(); }

    public int position() { return mVariant.position(); }

    public String toString() { return String.format("var(%s:%d %s>%s) tier(%s)",
            chromosome(), position(), mVariant.ref(), mVariant.alt(), mTier); }
}

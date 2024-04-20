package com.hartwig.hmftools.sage.candidate;

import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantTier;

public class Candidate
{
    private final VariantTier mTier;
    private VariantReadContext mReadContext;

    private int mMinNumberOfEvents;
    private int mReadContextSupport;

    public Candidate(final VariantTier tier, final VariantReadContext readContext, int minNumberOfEvents, int readContextSupport)
    {
        mTier = tier;
        mReadContext = readContext;
        mMinNumberOfEvents = minNumberOfEvents;
        mReadContextSupport = readContextSupport;
    }

    public static Candidate fromAltContext(final VariantTier tier, final AltContext altContext)
    {
        return new Candidate(tier, altContext.readContext(), altContext.minNumberOfEvents(), altContext.readContextSupport());
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
    public SimpleVariant variant() { return mReadContext.variant(); }
    public VariantReadContext readContext() { return mReadContext; }
    public int minNumberOfEvents() { return mMinNumberOfEvents; }
    public String chromosome() { return mReadContext.variant().Chromosome; }
    public int position() { return mReadContext.variant().Position; }

    public String toString() { return String.format("var(%s) tier(%s)", mReadContext.variant(), mTier); }
}

package com.hartwig.hmftools.sage.candidate;

import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.common.variant.VariantTier;

public class Candidate
{
    private final VariantTier mTier;
    private VariantReadContext mReadContext;

    private int mMinNumberOfEvents;
    private int mFullMatchSupport;
    private int mCoreMatchSupport;
    private int mLowQualInCoreCount;

    public Candidate(
            final VariantTier tier, final VariantReadContext readContext, int minNumberOfEvents, int fullMatchSupport)
    {
        mTier = tier;
        mReadContext = readContext;
        mMinNumberOfEvents = minNumberOfEvents;
        mFullMatchSupport = fullMatchSupport;
        mCoreMatchSupport = 0;
        mLowQualInCoreCount = 0;
    }

    public static Candidate fromAltContext(final VariantTier tier, final AltContext altContext)
    {
        Candidate candidate = new Candidate(
                tier, altContext.readContext(), altContext.minNumberOfEvents(), altContext.fullMatchSupport());

        candidate.setCoreMatchSupport(altContext.coreMatchSupport());
        candidate.setLowQualInCoreCount(altContext.lowQualInCoreCount());
        return candidate;
    }

    public void update(final AltContext altContext)
    {
        int altContextSupport = altContext.fullMatchSupport();
        if(altContextSupport > mFullMatchSupport)
        {
            mFullMatchSupport = altContextSupport;
            mReadContext = altContext.readContext();
            mMinNumberOfEvents = Math.min(mMinNumberOfEvents, altContext.minNumberOfEvents());
        }
    }

    public VariantTier tier() { return mTier; }
    public SimpleVariant variant() { return mReadContext.variant(); }
    public VariantReadContext readContext() { return mReadContext; }
    public int minNumberOfEvents() { return mMinNumberOfEvents; }
    public int fullMatchSupport() { return mFullMatchSupport; }
    public String chromosome() { return mReadContext.variant().Chromosome; }
    public int position() { return mReadContext.variant().Position; }

    public void setCoreMatchSupport(int count) { mCoreMatchSupport = count; }
    public int coreMatchSupport() { return mCoreMatchSupport; }

    public void setLowQualInCoreCount(int count) { mLowQualInCoreCount = count; }
    public int lowQualInCoreCount() { return mLowQualInCoreCount; }

    public String toString() { return String.format("var(%s) tier(%s)", mReadContext.variant(), mTier); }
}

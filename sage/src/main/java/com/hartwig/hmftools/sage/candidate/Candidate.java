package com.hartwig.hmftools.sage.candidate;

import com.google.common.annotations.VisibleForTesting;
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

    public Candidate(
            final VariantTier tier, final VariantReadContext readContext, int minNumberOfEvents, int fullMatchSupport,
            int coreMatchSupport)
    {
        mTier = tier;
        mReadContext = readContext;
        mMinNumberOfEvents = minNumberOfEvents;
        mFullMatchSupport = fullMatchSupport;
        mCoreMatchSupport = coreMatchSupport;
    }

    public Candidate(final VariantTier tier, final VariantReadContext readContext)
    {
        this(tier, readContext, 0, 0, 0);
    }

    public static Candidate fromAltCandidate(final VariantTier tier, final ReadContextCandidate readContextCandidate)
    {
        Candidate candidate = new Candidate(
                tier, readContextCandidate.readContext(), readContextCandidate.MinNumberOfEvents, readContextCandidate.FullMatch,
                readContextCandidate.CoreMatch);

        return candidate;
    }

    public void updateAltCandidate(final ReadContextCandidate altCandidate)
    {
        int altContextSupport = altCandidate.FullMatch;
        if(altContextSupport > mFullMatchSupport)
        {
            mFullMatchSupport = altContextSupport;
            mReadContext = altCandidate.readContext();
            mMinNumberOfEvents = Math.min(mMinNumberOfEvents, altCandidate.MinNumberOfEvents);
        }
    }

    public VariantTier tier() { return mTier; }
    public SimpleVariant variant() { return mReadContext.variant(); }
    public VariantReadContext readContext() { return mReadContext; }
    public String chromosome() { return mReadContext.variant().Chromosome; }
    public int position() { return mReadContext.variant().Position; }

    public int minNumberOfEvents() { return mMinNumberOfEvents; }
    public int fullMatchSupport() { return mFullMatchSupport; }
    public int coreMatchSupport() { return mCoreMatchSupport; }

    public String toString() { return String.format("var(%s) tier(%s)", mReadContext.variant(), mTier); }

    @VisibleForTesting
    public Candidate(final VariantTier tier, final VariantReadContext readContext, int minNumberOfEvents, int fullMatchSupport)
    {
        this(tier, readContext, minNumberOfEvents, fullMatchSupport, 0);
    }
}

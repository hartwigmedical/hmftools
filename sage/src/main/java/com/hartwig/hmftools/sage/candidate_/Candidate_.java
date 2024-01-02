package com.hartwig.hmftools.sage.candidate_;

import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantTier;

/**
 * Essentially a wrapper of AltContext + VarientTier.
 * Allows other AltContexts to be merged in.
 */
public class Candidate_
{
    private final VariantTier mTier;
    private final SimpleVariant mVariant;

    private int mMinNumberOfEvents;
    private int mReadContextSupport;
    private ReadContext_ mReadContext;

    public Candidate_(
            final VariantTier tier, final SimpleVariant variant, final ReadContext_ readContext, int minNumberOfEvents, int readContextSupport)
    {
        mTier = tier;
        mVariant = variant;
        mReadContext = readContext;
        mMinNumberOfEvents = minNumberOfEvents;
        mReadContextSupport = readContextSupport;
    }

    public static Candidate_ fromAltContext(final VariantTier tier, final AltContext_ altContext)
    {
        SimpleVariant variant = new SimpleVariant(altContext.chromosome(), altContext.position(), altContext.Ref, altContext.Alt);

        return new Candidate_(
                tier, variant, altContext.readContext(), altContext.minNumberOfEvents(), altContext.readContextSupport());
    }

    /**
     * Merges in another AltContext.
     * <p>
     * If the number of FULL + PARTIAL matches of altContext is larger than mReadContextSupport, then swap out these two, and the
     * readContext from altContext with mReadContext, and mMinNumberEvents is the min of its current value and
     * altContext.minNumberOfEvents().
     */
    public void update(final AltContext_ altContext)
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

    public SimpleVariant variant() { return mVariant; }

    public ReadContext_ readContext() { return mReadContext; }

    public int minNumberOfEvents() { return mMinNumberOfEvents; }

    public String chromosome() { return mVariant.chromosome(); }

    public int position() { return mVariant.position(); }

    public String toString() { return String.format("var(%s) tier(%s)", mVariant, mTier); }
}

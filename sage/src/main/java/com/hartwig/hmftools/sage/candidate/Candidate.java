package com.hartwig.hmftools.sage.candidate;

import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.common.VariantTier;

public class Candidate
{
    private final VariantTier mTier;
    private final VariantHotspot mVariant;

    private int mMaxDepth;
    private int mMinNumberOfEvents;
    private int mReadContextSupport;
    private ReadContext mReadContext;
    private int mRawSupportAlt;
    private int mRawBaseQualityAlt;

    public Candidate(
            final VariantTier tier, final VariantHotspot variant, final ReadContext readContext,
            int maxDepth, int minNumberOfEvents, int readContextSupport, int rawSupportAlt, int rawBaseQualAlt)
    {
        mTier = tier;
        mVariant = variant;
        mReadContext = readContext;
        mMaxDepth = maxDepth;
        mMinNumberOfEvents = minNumberOfEvents;
        mReadContextSupport = readContextSupport;
        mRawSupportAlt = rawSupportAlt;
        mRawBaseQualityAlt = rawBaseQualAlt;
    }

    public static Candidate fromAltContext(final VariantTier tier, final AltContext altContext)
    {
        return new Candidate(
                tier, ImmutableVariantHotspotImpl.builder().from(altContext).build(), altContext.readContext(),
                altContext.rawDepth(), altContext.minNumberOfEvents(), altContext.readContextSupport(),
                altContext.rawAltSupport(), altContext.rawAltBaseQuality());
    }

    public void update(final AltContext altContext)
    {
        int altContextSupport = altContext.readContextSupport();
        if(altContextSupport > mReadContextSupport)
        {
            mReadContextSupport = altContextSupport;
            mReadContext = altContext.readContext();
            mMinNumberOfEvents = Math.min(mMinNumberOfEvents, altContext.minNumberOfEvents());
            mRawSupportAlt = altContext.rawAltSupport();
            mRawBaseQualityAlt = altContext.rawAltBaseQuality();

        }

        mMaxDepth = Math.max(mMaxDepth, altContext.rawDepth());
    }

    public VariantTier tier() { return mTier; }
    public VariantHotspot variant() { return mVariant; }

    public ReadContext readContext() { return mReadContext; }

    public int maxReadDepth() { return mMaxDepth; }
    public int minNumberOfEvents() { return mMinNumberOfEvents; }
    public int rawSupportAlt() { return mRawSupportAlt; }
    public int rawBaseQualityAlt() { return mRawBaseQualityAlt; }

    public String chromosome() { return mVariant.chromosome(); }

    public int position() { return mVariant.position(); }

    public String toString() { return String.format("var(%s:%d %s>%s) tier(%s)",
            chromosome(), position(), mVariant.ref(), mVariant.alt(), mTier); }
}

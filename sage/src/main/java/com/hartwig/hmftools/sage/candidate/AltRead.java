package com.hartwig.hmftools.sage.candidate;

import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;

import htsjdk.samtools.SAMRecord;

public class AltRead
{
    public final String Ref;
    public final String Alt;
    public final int NumberOfEvents;
    public final boolean SufficientMapQuality;

    private final RefContext mRefContext;

    private final int mVariantReadIndex;
    private final SAMRecord mRead;

    public AltRead(
            final RefContext refContext, final String ref, final String alt, final int numberOfEvents,
            final boolean sufficientMapQuality, final SAMRecord read, final int variantReadIndex)
    {
        mRefContext = refContext;
        Ref = ref;
        Alt = alt;
        NumberOfEvents = numberOfEvents;
        SufficientMapQuality = sufficientMapQuality;

        mRead = read;
        mVariantReadIndex = variantReadIndex;
    }

    public int position()
    {
        return mRefContext.position();
    }

    public boolean isIndel() { return Ref.length() != Alt.length(); }

    public int length()
    {
        return Math.abs(Ref.length() - Alt.length());
    }

    public void updateRefContext(final VariantReadContextBuilder readContextBuilder, final RefSequence refSequence)
    {
        mRefContext.processAltRead(Ref, Alt, NumberOfEvents, mRead, mVariantReadIndex, readContextBuilder, refSequence);
    }

    public String toString() { return String.format("%d: %s>%s", mRefContext.Position, Ref, Alt); }
}

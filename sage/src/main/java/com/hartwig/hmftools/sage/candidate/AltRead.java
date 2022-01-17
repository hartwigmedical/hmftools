package com.hartwig.hmftools.sage.candidate;

import com.hartwig.hmftools.sage.common.ReadContext;

import org.jetbrains.annotations.Nullable;

public class AltRead
{
    public final String Ref;
    public final String Alt;
    public final int BaseQuality;
    public final int NumberOfEvents;
    public final boolean SufficientMapQuality;

    private final RefContext mRefContext;

    @Nullable
    private ReadContext mReadContext;

    public AltRead(
            final RefContext refContext, final String ref, final String alt, final int baseQuality, final int numberOfEvents,
            final boolean sufficientMapQuality, final ReadContext readContext)
    {
        mRefContext = refContext;
        Ref = ref;
        Alt = alt;
        BaseQuality = baseQuality;
        NumberOfEvents = numberOfEvents;
        SufficientMapQuality = sufficientMapQuality;

        mReadContext = readContext;
    }

    public boolean containsReadContext()
    {
        return mReadContext != null;
    }

    public int position()
    {
        return mRefContext.position();
    }

    public boolean isIndel()
    {
        return Ref.length() != Alt.length();
    }

    public int length()
    {
        return Math.abs(Ref.length() - Alt.length());
    }

    public int rightCoreIndex()
    {
        return mReadContext.readBasesRightCentreIndex();
    }
    public int leftCoreIndex()
    {
        return mReadContext.readBasesLeftCentreIndex();
    }

    public void extend(final AltRead other)
    {
        int leftIndex = Math.min(mReadContext.readBasesLeftCentreIndex(), other.mReadContext.readBasesLeftCentreIndex());
        int rightIndex = Math.max(mReadContext.readBasesRightCentreIndex(), other.mReadContext.readBasesRightCentreIndex());

        mReadContext.extendCore(leftIndex, rightIndex);
    }

    public void updateRefContext()
    {
        mRefContext.processAltRead(Ref, Alt, BaseQuality, SufficientMapQuality, NumberOfEvents, mReadContext);
    }

    public String toString() { return String.format("%s>%s", Ref, Alt); }
}

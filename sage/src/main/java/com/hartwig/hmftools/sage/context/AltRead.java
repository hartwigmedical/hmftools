package com.hartwig.hmftools.sage.context;

import com.hartwig.hmftools.sage.read.ReadContext;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class AltRead
{
    @Nullable
    public final RefContext RefContext;
    public final String Ref;
    public final String Alt;
    public final int BaseQuality;
    public final int NumberOfEvents;
    public final boolean SufficientMapQuality;
    
    private ReadContext mReadContext;

    public AltRead(
            final RefContext refContext, final String ref, final String alt, final int baseQuality, final int numberOfEvents,
            final boolean sufficientMapQuality, final ReadContext readContext)
    {
        RefContext = refContext;
        Ref = ref;
        Alt = alt;
        BaseQuality = baseQuality;
        mReadContext = readContext;
        NumberOfEvents = numberOfEvents;
        SufficientMapQuality = sufficientMapQuality;
    }

    public boolean containsReadContext()
    {
        return mReadContext != null;
    }

    public int position()
    {
        return RefContext.position();
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

    public void extend(@NotNull final AltRead other)
    {
        assert (mReadContext != null);
        assert (other.mReadContext != null);

        int leftIndex = Math.min(mReadContext.readBasesLeftCentreIndex(), other.mReadContext.readBasesLeftCentreIndex());
        int rightIndex = Math.max(mReadContext.readBasesRightCentreIndex(), other.mReadContext.readBasesRightCentreIndex());

        mReadContext = mReadContext.extend(leftIndex, rightIndex);
    }

    public void updateRefContext()
    {
        RefContext.altRead(Ref, Alt, BaseQuality, SufficientMapQuality, NumberOfEvents, mReadContext);
    }

}

package com.hartwig.hmftools.sage.candidate_;

import org.jetbrains.annotations.Nullable;

/**
 * Contains a ReadContext and variant info for an alt read.
 * Also links back to a ref context from which you can update with this AltRead.
 * <p>
 * Contains some other counters.
 */
public class AltRead_
{
    public final String Ref;
    public final String Alt;
    public final int BaseQuality;
    public final int NumberOfEvents;
    public final boolean SufficientMapQuality;

    private final RefContext_ mRefContext;

    @Nullable
    private ReadContext_ mReadContext;

    public AltRead_(
            final RefContext_ refContext, final String ref, final String alt, final int baseQuality, final int numberOfEvents,
            final boolean sufficientMapQuality, @Nullable final ReadContext_ readContext)
    {
        mRefContext = refContext;
        Ref = ref;
        Alt = alt;
        BaseQuality = baseQuality;
        NumberOfEvents = numberOfEvents;
        SufficientMapQuality = sufficientMapQuality;

        mReadContext = readContext;
    }

    /**
     * Triggers an update to the ref context.
     */
    public void updateRefContext()
    {
        mRefContext.processAltRead(Ref, Alt, BaseQuality, NumberOfEvents, mReadContext);
    }

    /**
     * Extend the core of the read context based on another AltRead.
     * @param other
     */
    public void extend(final AltRead_ other)
    {
        int leftIndex = Math.min(mReadContext.readBasesLeftCentreIndex(), other.mReadContext.readBasesLeftCentreIndex());
        int rightIndex = Math.max(mReadContext.readBasesRightCentreIndex(), other.mReadContext.readBasesRightCentreIndex());

        mReadContext.extendCore(leftIndex, rightIndex);
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

    /**
     * Indel length.
     */
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

    public String toString() { return String.format("%d: %s>%s", position(), Ref, Alt); }
}

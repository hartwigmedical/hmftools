package com.hartwig.hmftools.sage.context;

import com.hartwig.hmftools.sage.read.ReadContext;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class AltRead
{

    @Nullable
    private ReadContext readContext;
    private final RefContext refContext;
    private final String ref;
    private final String alt;
    private final int baseQuality;
    private final int numberOfEvents;
    private final boolean sufficientMapQuality;

    AltRead(final RefContext refContext, final String ref, final String alt, final int baseQuality, final int numberOfEvents,
            final boolean sufficientMapQuality,
            @Nullable final ReadContext readContext)
    {
        this.refContext = refContext;
        this.ref = ref;
        this.alt = alt;
        this.baseQuality = baseQuality;
        this.readContext = readContext;
        this.numberOfEvents = numberOfEvents;
        this.sufficientMapQuality = sufficientMapQuality;
    }

    public boolean containsReadContext()
    {
        return readContext != null;
    }

    public long position()
    {
        return refContext.position();
    }

    public boolean isIndel()
    {
        return ref.length() != alt.length();
    }

    public int length()
    {
        return Math.abs(ref.length() - alt.length());
    }

    public int rightCoreIndex()
    {
        return readContext.readBasesRightCentreIndex();
    }

    public int leftCoreIndex()
    {
        return readContext.readBasesLeftCentreIndex();
    }

    public void extend(@NotNull final AltRead other)
    {
        assert (readContext != null);
        assert (other.readContext != null);

        int leftIndex = Math.min(readContext.readBasesLeftCentreIndex(), other.readContext.readBasesLeftCentreIndex());
        int rightIndex = Math.max(readContext.readBasesRightCentreIndex(), other.readContext.readBasesRightCentreIndex());

        readContext = readContext.extend(leftIndex, rightIndex);
    }

    public void updateRefContext()
    {
        refContext.altRead(ref, alt, baseQuality, sufficientMapQuality, numberOfEvents, readContext);
    }

}

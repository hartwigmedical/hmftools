package com.hartwig.hmftools.sage.read;

import org.jetbrains.annotations.NotNull;

public class RawContext
{

    private final int readIndex;
    private final boolean readIndexInSoftClip;
    private final boolean readIndexInDelete;
    private final boolean readIndexInSkipped;
    private final boolean indelAtPosition;
    private final boolean altSupport;
    private final boolean refSupport;
    private final boolean depthSupport;
    private final int altQuality;
    private final int refQuality;

    @NotNull
    static RawContext inSoftClip(final int readIndex)
    {
        return new RawContext(readIndex, false, false, true, false, false, false, false, 0, 0);
    }

    @NotNull
    static RawContext inDelete(final int readIndex)
    {
        return new RawContext(readIndex, true, false, false, false, false, false, false, 0, 0);
    }

    @NotNull
    static RawContext inSkipped(final int readIndex)
    {
        return new RawContext(readIndex, false, true, false, false, false, false, false, 0, 0);
    }

    @NotNull
    static RawContext indel(final int readIndex, final boolean altSupport, final int quality)
    {
        return new RawContext(readIndex, false, false, false, true, altSupport, false, true, altSupport ? quality : 0, 0);
    }

    @NotNull
    static RawContext alignment(final int readIndex, final boolean altSupport, final boolean refSupport, final int quality)
    {
        return new RawContext(readIndex,
                false,
                false,
                false,
                false,
                altSupport,
                refSupport,
                true,
                altSupport ? quality : 0,
                refSupport ? quality : 0);
    }

    public RawContext(final int readIndex, final boolean readIndexInDelete, final boolean readIndexInSkipped,
            final boolean readIndexInSoftClip, final boolean indelAtPosition, final boolean altSupport, final boolean refSupport,
            final boolean depthSupport, final int altQuality, final int refQuality)
    {
        this.readIndex = readIndex;
        this.readIndexInDelete = readIndexInDelete;
        this.readIndexInSkipped = readIndexInSkipped;
        this.readIndexInSoftClip = readIndexInSoftClip;
        this.altSupport = altSupport;
        this.refSupport = refSupport;
        this.depthSupport = depthSupport;
        this.altQuality = altQuality;
        this.refQuality = refQuality;
        this.indelAtPosition = indelAtPosition;
    }

    public int readIndex()
    {
        return readIndex;
    }

    public boolean isIndelAtPosition()
    {
        return indelAtPosition;
    }

    public boolean isReadIndexInDelete()
    {
        return readIndexInDelete;
    }

    public boolean isReadIndexInSoftClip()
    {
        return readIndexInSoftClip;
    }

    public boolean isReadIndexInSkipped()
    {
        return readIndexInSkipped;
    }

    public boolean isAltSupport()
    {
        return altSupport;
    }

    public boolean isRefSupport()
    {
        return refSupport;
    }

    public boolean isDepthSupport()
    {
        return depthSupport;
    }

    public int altQuality()
    {
        return altQuality;
    }

    public int refQuality()
    {
        return refQuality;
    }
}

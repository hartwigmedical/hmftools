package com.hartwig.hmftools.sage.read;

public class RawContext
{
    public final int ReadIndex;
    public final boolean ReadIndexInSoftClip;
    public final boolean ReadIndexInDelete;
    public final boolean ReadIndexInSkipped;
    public final boolean IndelAtPosition;
    public final boolean AltSupport;
    public final boolean RefSupport;
    public final boolean DepthSupport;
    public final int AltQuality;
    public final int RefQuality;

    static RawContext inSoftClip(final int readIndex)
    {
        return new RawContext(
                readIndex, false, false, true, false,
                false, false, false, 0, 0);
    }

    static RawContext inDelete(final int readIndex)
    {
        return new RawContext(
                readIndex, true, false, false, false,
                false, false, false, 0, 0);
    }

    static RawContext inSkipped(final int readIndex)
    {
        return new RawContext(
                readIndex, false, true, false, false,
                false, false, false, 0, 0);
    }

    static RawContext indel(final int readIndex, final boolean altSupport, final int quality)
    {
        return new RawContext(
                readIndex, false, false, false, true,
                altSupport, false, true, altSupport ? quality : 0, 0);
    }

    static RawContext alignment(final int readIndex, final boolean altSupport, final boolean refSupport, final int quality)
    {
        return new RawContext(
                readIndex, false, false, false, false,
                altSupport, refSupport, true, altSupport ? quality : 0, refSupport ? quality : 0);
    }

    public RawContext(
            final int readIndex, final boolean readIndexInDelete, final boolean readIndexInSkipped,
            final boolean readIndexInSoftClip, final boolean indelAtPosition, final boolean altSupport, final boolean refSupport,
            final boolean depthSupport, final int altQuality, final int refQuality)
    {
        ReadIndex = readIndex;
        ReadIndexInDelete = readIndexInDelete;
        ReadIndexInSkipped = readIndexInSkipped;
        ReadIndexInSoftClip = readIndexInSoftClip;
        AltSupport = altSupport;
        RefSupport = refSupport;
        DepthSupport = depthSupport;
        AltQuality = altQuality;
        RefQuality = refQuality;
        IndelAtPosition = indelAtPosition;
    }
}

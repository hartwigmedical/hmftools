package com.hartwig.hmftools.sage.evidence;

import com.hartwig.hmftools.common.samtools.CigarTraversal;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import htsjdk.samtools.SAMRecord;

public class RawContext
{
    public final int ReadIndex;
    public final boolean ReadIndexInSoftClip;
    public final boolean ReadIndexInDelete;
    public final boolean ReadIndexInSkipped;
    public final boolean AltSupport;
    public final boolean RefSupport;
    public final boolean DepthSupport;
    public final int AltQuality;
    public final int RefQuality;

    protected static final RawContext INVALID_CONTEXT = new RawContext(
            -1, false, false, false,
            false, false, false, 0, 0);

    public RawContext(
            final int readIndex, final boolean readIndexInDelete, final boolean readIndexInSkipped,
            final boolean readIndexInSoftClip, final boolean altSupport, final boolean refSupport,
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
    }

    public static RawContext create(final VariantHotspot variant, final SAMRecord record)
    {
        RawContextCigarHandler handler = new RawContextCigarHandler(variant);
        CigarTraversal.traverseCigar(record, handler);
        RawContext result = handler.result();
        return result == null ? INVALID_CONTEXT : result;
    }

    static RawContext inSoftClip(final int readIndex, final boolean altSupport, final int quality)
    {
        return new RawContext(
                readIndex, false, false, true,
                altSupport, false, false, altSupport ? quality : 0, 0);
    }

    static RawContext inDelete(final int readIndex)
    {
        return new RawContext(
                readIndex, true, false, false,
                false, false, false, 0, 0);
    }

    static RawContext inSkipped(final int readIndex)
    {
        return new RawContext(
                readIndex, false, true, false,
                false, false, false, 0, 0);
    }

    static RawContext indel(final int readIndex, final boolean altSupport, final int quality)
    {
        return new RawContext(
                readIndex, false, false, false,
                altSupport, false, true, altSupport ? quality : 0, 0);
    }

    static RawContext alignment(final int readIndex, final boolean altSupport, final boolean refSupport, final int quality)
    {
        return new RawContext(
                readIndex, false, false, false,
                altSupport, refSupport, true, altSupport ? quality : 0, refSupport ? quality : 0);
    }
}

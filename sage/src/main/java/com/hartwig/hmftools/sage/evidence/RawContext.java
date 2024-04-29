package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

import com.hartwig.hmftools.common.bam.CigarHandler;
import com.hartwig.hmftools.sage.common.SimpleVariant;

import htsjdk.samtools.SAMRecord;

public class RawContext
{
    public final int ReadIndex;
    public final boolean ReadIndexInSoftClip;
    public final boolean ReadIndexInDelete;
    public final boolean ReadIndexInSkipped;
    public final int BaseQuality;

    protected static final RawContext INVALID_CONTEXT = new RawContext(
            -1, false, false, false,0);

    public RawContext(
            final int readIndex, final boolean readIndexInDelete, final boolean readIndexInSkipped,
            final boolean readIndexInSoftClip, final int baseQuality)
    {
        ReadIndex = readIndex;
        ReadIndexInDelete = readIndexInDelete;
        ReadIndexInSkipped = readIndexInSkipped;
        ReadIndexInSoftClip = readIndexInSoftClip;
        BaseQuality = baseQuality;
    }

    public static RawContext create(final SimpleVariant variant, final SAMRecord record)
    {
        RawContextCigarHandler handler = new RawContextCigarHandler(variant);
        CigarHandler.traverseCigar(record, handler);
        RawContext result = handler.result();
        return result == null ? INVALID_CONTEXT : result;
    }

    public String toString()
    {
        return format("index(%d) sc(%s) del(%s) skip(%s) bq(%d)",
                ReadIndex, ReadIndexInSoftClip, ReadIndexInDelete, ReadIndexInSkipped, BaseQuality);
    }
}

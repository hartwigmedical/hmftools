package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

import com.hartwig.hmftools.common.bam.CigarHandler;
import com.hartwig.hmftools.sage.common.SimpleVariant;

import htsjdk.samtools.SAMRecord;

public class RawContext
{
    public final int ReadVariantIndex;
    public final VariantReadPositionType PositionType;

    protected static final RawContext INVALID_CONTEXT = new RawContext(-1, VariantReadPositionType.NONE);

    public RawContext(
            final int readVariantIndex, final VariantReadPositionType positionType)
    {
        ReadVariantIndex = readVariantIndex;
        PositionType = positionType;
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
        return format("readVarIndex(%d) posType(%s)", ReadVariantIndex, PositionType);
    }
}

package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.common.ReadCigarInfo.getReadIndexFromPosition;

import com.hartwig.hmftools.sage.common.VariantReadContext;

import htsjdk.samtools.SAMRecord;

public class Realignment
{
    public static final int INVALID_INDEX = -1;

    public static int realignedReadIndexPosition(final VariantReadContext readContext, final SAMRecord record)
    {
        int variantCoreEndPosition = readContext.CorePositionEnd;

        if(variantCoreEndPosition < record.getAlignmentStart() || variantCoreEndPosition > record.getAlignmentEnd())
            return INVALID_INDEX;

        int coreEndReadIndex = getReadIndexFromPosition(record.getAlignmentStart(), record.getCigar().getCigarElements(), variantCoreEndPosition);

        // convert back to the variant's index location
        int adjustedReadIndex = coreEndReadIndex - readContext.rightCoreLength();
        return adjustedReadIndex;
    }
}

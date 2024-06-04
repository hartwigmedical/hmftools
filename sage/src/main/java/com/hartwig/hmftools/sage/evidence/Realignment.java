package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.common.ReadCigarInfo.getReadIndexFromPosition;
import static com.hartwig.hmftools.sage.evidence.JitterMatch.checkJitter;
import static com.hartwig.hmftools.sage.evidence.RealignedType.EXACT;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.common.ReadContextMatcher;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import htsjdk.samtools.SAMRecord;

public class Realignment
{
    public static RealignedType checkRealignment(
            final VariantReadContext readContext, final ReadContextMatcher readContextMatcher, final SAMRecord record, final int readIndex)
    {
        // the read index corresponding to the ref position at the end of the core
        int realignedReadIndex = Realignment.realignedReadIndexPosition(readContext, record);

        if(readIndex == realignedReadIndex)
            return RealignedType.NONE;

        if(realignedReadIndex < 0 || realignedReadIndex >= record.getReadBases().length)
            return RealignedType.NONE;

        ReadContextMatch match = readContextMatcher.determineReadMatch(record, realignedReadIndex);

        if(match == ReadContextMatch.FULL || match == ReadContextMatch.PARTIAL_CORE)
            return RealignedType.EXACT;

        // otherwise check jitter
        JitterMatch jitterMatch = checkJitter(readContext, record, realignedReadIndex);

        if(jitterMatch == JitterMatch.SHORTENED)
            return RealignedType.SHORTENED;
        else if(jitterMatch == JitterMatch.LENGTHENED)
            return RealignedType.SHORTENED;

        return RealignedType.NONE;
    }

    public static final int INVALID_INDEX = -1;

    @VisibleForTesting
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

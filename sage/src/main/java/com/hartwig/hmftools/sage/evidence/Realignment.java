package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.sage.common.ReadCigarInfo.getReadIndexFromPosition;
import static com.hartwig.hmftools.sage.evidence.JitterMatch.checkJitter;
import static com.hartwig.hmftools.sage.evidence.RealignedType.EXACT;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.common.ReadContextMatcher;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import htsjdk.samtools.SAMRecord;

public class Realignment
{
    public static RealignedType checkRealignment(
            final VariantReadContext readContext, final ReadContextMatcher readContextMatcher, final SAMRecord record,
            final int readIndex, final int realignedReadIndex)
    {
        // the read index corresponding to the ref position at the end of the core
        if(readIndex == realignedReadIndex)
            return RealignedType.NONE;

        if(realignedReadIndex < 0 || realignedReadIndex >= record.getReadBases().length)
            return RealignedType.NONE;

        ReadContextMatch match = readContextMatcher.determineReadMatch(
                record.getReadBases(), record.getBaseQualities(), realignedReadIndex, true);

        if(match == ReadContextMatch.FULL || match == ReadContextMatch.PARTIAL_CORE)
            return RealignedType.EXACT;

        // otherwise check jitter
        JitterMatch jitterMatch = checkJitter(readContext, record, realignedReadIndex);

        if(jitterMatch == JitterMatch.SHORTENED)
            return RealignedType.SHORTENED;
        else if(jitterMatch == JitterMatch.LENGTHENED)
            return RealignedType.LENGTHENED;

        return RealignedType.NONE;
    }

    public static boolean considerRealignedDel(final SimpleVariant variant, final SAMRecord record)
    {
        if(!variant.isDelete())
            return false;

        int posDiff = record.getAlignmentStart() - variant.Position;
        return posDiff <= variant.refLength();
    }

    public static final int INVALID_INDEX = -1;

    public static int realignedReadIndexPosition(final VariantReadContext readContext, final SAMRecord record)
    {
        int variantCoreEndPosition = readContext.CorePositionEnd;

        if(variantCoreEndPosition < record.getAlignmentStart() || variantCoreEndPosition > record.getAlignmentEnd())
            return INVALID_INDEX;

        int coreEndReadIndex;
        if(variantCoreEndPosition <= record.getAlignmentStart())
        {
            int leftClipLength = CigarUtils.leftSoftClipLength(record);
            coreEndReadIndex = leftClipLength - (record.getAlignmentStart() - variantCoreEndPosition);
        }
        else
        {
            coreEndReadIndex = getReadIndexFromPosition(
                    record.getAlignmentStart(), record.getCigar().getCigarElements(), variantCoreEndPosition);
        }

        // convert back to the variant's index location
        int adjustedReadIndex = coreEndReadIndex - readContext.rightCoreLength();
        return adjustedReadIndex;
    }
}

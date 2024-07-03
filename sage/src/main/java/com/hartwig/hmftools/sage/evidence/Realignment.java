package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.bam.CigarUtils.getReadIndexFromPosition;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.INVALID_READ_INDEX;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.evidence.JitterMatch.checkJitter;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.common.ReadContextMatcher;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import htsjdk.samtools.SAMRecord;

public class Realignment
{
    public static RealignedType checkRealignment(
            final VariantReadContext readContext, final ReadContextMatcher readContextMatcher, final SAMRecord record,
            final int readIndex, int realignedReadIndex, final SplitReadSegment splitReadSegment)
    {
        // the read index corresponding to the ref position at the end of the core
        if(readIndex == realignedReadIndex)
            return RealignedType.NONE;

        ReadContextMatch match = NONE;

        int realignmentOffset = realignedReadIndex - readIndex;

        if(splitReadSegment != null)
        {
            realignedReadIndex -= splitReadSegment.SegmentIndexStart;

            if(realignedReadIndex < 0 || realignedReadIndex >= splitReadSegment.length())
                return RealignedType.NONE;

            match = checkMatch(
                    readContextMatcher, splitReadSegment.ReadBases, splitReadSegment.ReadQuals, realignedReadIndex, realignmentOffset);
        }
        else
        {
            if(realignedReadIndex < 0 || realignedReadIndex >= record.getReadBases().length)
                return RealignedType.NONE;

            match = checkMatch(
                    readContextMatcher, record.getReadBases(), record.getBaseQualities(), realignedReadIndex, realignmentOffset);
        }

        if(match == ReadContextMatch.FULL || match == ReadContextMatch.PARTIAL_CORE)
            return RealignedType.EXACT;

        // otherwise check jitter
        JitterMatch jitterMatch = checkJitter(readContext, readContextMatcher, record, realignedReadIndex);

        if(jitterMatch == JitterMatch.SHORTENED)
            return RealignedType.SHORTENED;
        else if(jitterMatch == JitterMatch.LENGTHENED)
            return RealignedType.LENGTHENED;

        return RealignedType.NONE;
    }

    private static ReadContextMatch checkMatch(
            final ReadContextMatcher readContextMatcher, final byte[] readBases, final byte[] readQuals, final int readVarIndex,
            int realignmentOffset)
    {
        readContextMatcher.setRealignmentIndexOffset(realignmentOffset);

        ReadContextMatch match = readContextMatcher.determineReadMatch(readBases, readQuals, readVarIndex, true);

        readContextMatcher.clearRealignmentIndexOffset();

        return match;
    }

    public static boolean considerRealignedDel(final SAMRecord record, final int minPositionVsRead)
    {
        // eg a DEL at position 100 with 4 bases deleted (ref length 5) needs to consider reads starting as late as 104
        int unclippedStartPos = record.getAlignmentStart() - CigarUtils.leftSoftClipLength(record);
        return minPositionVsRead >= unclippedStartPos;
    }

    public static final int INVALID_INDEX = -1;

    public static int realignedReadIndexPosition(final VariantReadContext readContext, final SAMRecord record)
    {
        int variantCoreEndPosition = readContext.CorePositionEnd;

        int coreEndReadIndex = getReadIndexFromPosition(
                record.getAlignmentStart(), record.getCigar().getCigarElements(), variantCoreEndPosition,
                false, true);

        if(coreEndReadIndex == INVALID_READ_INDEX)
            return INVALID_INDEX;

        // convert back to the variant's index location
        int adjustedReadIndex = coreEndReadIndex - readContext.rightCoreLength();
        return adjustedReadIndex;
    }
}

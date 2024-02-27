package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;
import static com.hartwig.hmftools.sage.evidence.RealignedContext.NONE;
import static com.hartwig.hmftools.sage.evidence.RealignedType.EXACT;
import static com.hartwig.hmftools.sage.evidence.RealignedType.LENGTHENED;
import static com.hartwig.hmftools.sage.evidence.RealignedType.SHORTENED;

import com.hartwig.hmftools.common.variant.repeat.RepeatContextFactory;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.ReadContext;

import org.jetbrains.annotations.Nullable;

public class Realignment
{
    private static final int MIN_REPEAT_COUNT = 4;
    public static final int MAX_REPEAT_SIZE = 5;

    private static final Repeat NO_REPEAT = new Repeat(0, 0);

    public static RealignedContext realignedAroundIndex(
            final ReadContext readContext, final int otherBaseIndex, final byte[] otherBases, final byte[] otherBaseQuals, int maxSize)
    {
        int baseStartIndex = readContext.readBasesLeftFlankIndex();
        int baseEndIndex = readContext.readBasesRightFlankIndex();

        int leftOffset = readContext.readBasesPositionIndex() - baseStartIndex;
        int otherStartIndex = otherBaseIndex - leftOffset;

        return realigned(readContext, baseStartIndex, baseEndIndex, readContext.readBases(), otherStartIndex, otherBases, otherBaseQuals, maxSize);
    }

    public static RealignedContext realigned(
            final ReadContext readContext, int baseStartIndex, int baseEndIndex, final byte[] bases, final int otherBaseIndex,
            final byte[] otherBases, final byte[] otherBaseQuals, int maxDistance)
    {
        if(otherBaseIndex >= 0)
        {
            final RealignedContext context = realigned(
                    readContext, baseStartIndex, baseEndIndex, bases, otherBaseIndex, otherBases, otherBaseQuals);

            if(context.Type != RealignedType.NONE)
                return context;
        }

        RealignedContext result = NONE;
        for(int i = -maxDistance; i <= maxDistance; i++)
        {
            int otherBaseIndexWithOffset = otherBaseIndex + i;
            if(i != 0 && otherBaseIndexWithOffset >= 0)
            {
                final RealignedContext context = realigned(
                        readContext, baseStartIndex, baseEndIndex, bases, otherBaseIndexWithOffset, otherBases, otherBaseQuals);
                if(context.Type != RealignedType.NONE)
                {
                    if(context.Type == RealignedType.EXACT)
                        return context;

                    result = context;
                }
            }
        }

        return result;
    }

    public static RealignedContext realigned(final ReadContext readContext, int baseStartIndex, int baseEndIndex, final byte[] bases,
            int otherIndex, byte[] otherBases,
            final byte[] otherBaseQuals)
    {
        int exactLength = baseEndIndex - baseStartIndex + 1;

        int matchingBases = matchingBasesFromLeft(readContext, baseStartIndex, baseEndIndex, bases, otherIndex, otherBases, otherBaseQuals);

        if(matchingBases == exactLength)
        {
            return new RealignedContext(EXACT, matchingBases, otherIndex);
        }

        if(matchingBases < MIN_REPEAT_COUNT)
            return NONE;

        int baseNextIndex = baseStartIndex + matchingBases;
        int otherNextIndex = otherIndex + matchingBases;

        final Repeat repeat = repeatCount(otherNextIndex, otherBases);
        int repeatLength = repeat.RepeatLength;

        if(repeatLength == 0)
            return NONE;

        int matchingBasesShortened = matchingBasesFromLeft(
                readContext, baseNextIndex + repeatLength, baseEndIndex, bases, otherNextIndex, otherBases, otherBaseQuals);
        if(matchingBasesShortened > 0 && matchingBases + matchingBasesShortened == exactLength - repeatLength)
        {
            return new RealignedContext(
                    SHORTENED, matchingBasesShortened, otherNextIndex, repeat.RepeatCount, repeatLength, otherIndex, matchingBases);
        }

        int matchingBasesLengthened = matchingBasesFromLeft(
                readContext, baseNextIndex - repeatLength, baseEndIndex, bases, otherNextIndex, otherBases, otherBaseQuals);
        if(matchingBasesLengthened > 0 && matchingBases + matchingBasesLengthened == exactLength + repeatLength)
        {
            return new RealignedContext(
                    LENGTHENED, matchingBasesLengthened, otherNextIndex, repeat.RepeatCount + 1, repeatLength, otherIndex, matchingBases);
        }

        return NONE;
    }

    private static int matchingBasesFromLeft(@Nullable final ReadContext readContext, int startIndex, int endIndex, byte[] bases,
            int otherStartIndex, byte[] otherBases, byte[] otherBaseQuals)
    {
        if(startIndex < 0)
            return 0;

        int maxLength = Math.min(endIndex - startIndex + 1, otherBases.length - otherStartIndex);

        IndexedBases indexedBases = readContext != null ? readContext.indexedBases() : null;
        for(int i = 0; i < maxLength; i++)
        {
            if(bases[startIndex + i] != otherBases[otherStartIndex + i])
            {
                boolean withinCoreAndFlanks = indexedBases == null ? false : positionWithin(startIndex + i, indexedBases.LeftFlankIndex, indexedBases.RightFlankIndex);
                boolean withinCore = indexedBases == null ? false : positionWithin(startIndex + i, indexedBases.LeftCoreIndex, indexedBases.RightCoreIndex);
                boolean withinFlanks = withinCoreAndFlanks && !withinCore;
                if(withinFlanks && isLowBaseQual(otherBaseQuals, otherStartIndex + i))
                {
                    continue;
                }

                return i;
            }
        }

        return maxLength;
    }

    private static boolean isLowBaseQual(final byte[] baseQualities, int bqIndex)
    {
        if(baseQualities == null)
        {
            return false;
        }

        if(bqIndex < 0 || bqIndex >= baseQualities.length)
        {
            return false;
        }

        return baseQualities[bqIndex] < MATCHING_BASE_QUALITY;
    }

    private static Repeat repeatCount(int index, byte[] bases)
    {
        for(int i = 1; i <= MAX_REPEAT_SIZE; i++)
        {
            int repeats = RepeatContextFactory.backwardRepeats(index - i, i, bases) + 1;

            if(repeats >= MIN_REPEAT_COUNT)
                return new Repeat(i, repeats);
        }

        return NO_REPEAT;
    }

    private static class Repeat
    {
        public final int RepeatLength;
        public final int RepeatCount;

        Repeat(final int repeatLength, final int repeatCount)
        {
            RepeatLength = repeatLength;
            RepeatCount = repeatCount;
        }
    }
}

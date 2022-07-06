package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.evidence.RealignedContext.NONE;
import static com.hartwig.hmftools.sage.evidence.RealignedType.EXACT;
import static com.hartwig.hmftools.sage.evidence.RealignedType.LENGTHENED;
import static com.hartwig.hmftools.sage.evidence.RealignedType.SHORTENED;

import com.hartwig.hmftools.common.variant.repeat.RepeatContextFactory;
import com.hartwig.hmftools.sage.common.ReadContext;

public class Realigned
{
    private static final int MIN_REPEAT_COUNT = 4;
    public static final int MAX_REPEAT_SIZE = 5;

    private static final Repeat NO_REPEAT = new Repeat(0, 0);

    public static RealignedContext realignedAroundIndex(
            final ReadContext readContext, final int otherBaseIndex, final byte[] otherBases, int maxSize)
    {
        int baseStartIndex = readContext.readBasesLeftFlankIndex();
        int baseEndIndex = readContext.readBasesRightFlankIndex();

        int leftOffset = readContext.readBasesPositionIndex() - baseStartIndex;
        int otherStartIndex = otherBaseIndex - leftOffset;

        return realigned(baseStartIndex, baseEndIndex, readContext.readBases(), otherStartIndex, otherBases, maxSize);
    }

    public static RealignedContext realigned(
            int baseStartIndex, int baseEndIndex, final byte[] bases, final int otherBaseIndex, final byte[] otherBases, int maxDistance)
    {
        if(otherBaseIndex >= 0)
        {
            final RealignedContext context = realigned(baseStartIndex, baseEndIndex, bases, otherBaseIndex, otherBases);

            if(context.Type != RealignedType.NONE)
                return context;
        }

        RealignedContext result = NONE;
        for(int i = -maxDistance; i <= maxDistance; i++)
        {
            int otherBaseIndexWithOffset = otherBaseIndex + i;
            if(i != 0 && otherBaseIndexWithOffset >= 0)
            {
                final RealignedContext context = realigned(baseStartIndex, baseEndIndex, bases, otherBaseIndexWithOffset, otherBases);
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

    public static RealignedContext realigned(int baseStartIndex, int baseEndIndex, final byte[] bases, int otherIndex, byte[] otherBases)
    {
        int exactLength = baseEndIndex - baseStartIndex + 1;

        int matchingBases = matchingBasesFromLeft(baseStartIndex, baseEndIndex, bases, otherIndex, otherBases);

        if(matchingBases == exactLength)
        {
            return new RealignedContext(EXACT, matchingBases, otherIndex, 0, 0);
        }

        if(matchingBases < MIN_REPEAT_COUNT)
            return NONE;

        int baseNextIndex = baseStartIndex + matchingBases;
        int otherNextIndex = otherIndex + matchingBases;

        final Repeat repeat = repeatCount(otherNextIndex, otherBases);
        int repeatLength = repeat.RepeatLength;

        if(repeatLength == 0)
            return NONE;

        int matchingBasesShortened = matchingBasesFromLeft(baseNextIndex + repeatLength, baseEndIndex, bases, otherNextIndex, otherBases);
        if(matchingBasesShortened > 0 && matchingBases + matchingBasesShortened == exactLength - repeatLength)
        {
            return new RealignedContext(SHORTENED, matchingBasesShortened, otherNextIndex, repeat.RepeatCount, repeatLength);
        }

        int matchingBasesLengthened = matchingBasesFromLeft(baseNextIndex - repeatLength, baseEndIndex, bases, otherNextIndex, otherBases);
        if(matchingBasesLengthened > 0 && matchingBases + matchingBasesLengthened == exactLength + repeatLength)
        {
            return new RealignedContext(LENGTHENED, matchingBasesLengthened, otherNextIndex, repeat.RepeatCount + 1, repeatLength);
        }

        return NONE;
    }

    private static int matchingBasesFromLeft(int startIndex, int endIndex, byte[] bases, int otherStartIndex, byte[] otherBases)
    {
        if(startIndex < 0)
            return 0;

        int maxLength = Math.min(endIndex - startIndex + 1, otherBases.length - otherStartIndex);

        for(int i = 0; i < maxLength; i++)
        {
            if(bases[startIndex + i] != otherBases[otherStartIndex + i])
                return i;
        }

        return maxLength;
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

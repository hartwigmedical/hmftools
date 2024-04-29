package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.max;

import static com.hartwig.hmftools.sage.evidence.RealignedContext.NONE;
import static com.hartwig.hmftools.sage.evidence.RealignedType.EXACT;
import static com.hartwig.hmftools.sage.evidence.RealignedType.LENGTHENED;
import static com.hartwig.hmftools.sage.evidence.RealignedType.SHORTENED;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;

import com.hartwig.hmftools.sage.common.VariantReadContext;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class Realignment
{
    public static final int INVALID_INDEX = -1;

    public static int realignedReadIndexPosition(final VariantReadContext readContext, final SAMRecord record)
    {
        int variantCoreEndPosition = readContext.CorePositionEnd;

        if(variantCoreEndPosition < record.getAlignmentStart() || variantCoreEndPosition > record.getAlignmentEnd())
            return INVALID_INDEX;

        int refPosition = record.getAlignmentStart();
        int index = 0;

        for(CigarElement element : record.getCigar().getCigarElements())
        {
            if(refPosition + element.getLength() >= variantCoreEndPosition && element.getOperator().consumesReferenceBases())
            {
                if(element.getOperator().consumesReadBases())
                    index += variantCoreEndPosition - refPosition;

                break;
            }

            if(element.getOperator().consumesReferenceBases())
                refPosition += element.getLength();

            if(element.getOperator().consumesReadBases())
                index += element.getLength();
        }

        // convert back to the variant's index location
        int adjustedReadIndex = index - readContext.rightCoreLength();
        return adjustedReadIndex;
    }

    public static RealignedContext realignedAroundIndex(final VariantReadContext readContext, final int readIndex, final SAMRecord record)
    {
        int baseStartIndex = 0;
        int baseEndIndex = readContext.totalLength() - 1;

        int leftOffset = readContext.leftLength();
        int otherStartIndex = readIndex - leftOffset;

        int maxAlignDistance = getMaxRealignDistance(readContext, record);

        return realigned(baseStartIndex, baseEndIndex, readContext.ReadBases, otherStartIndex, record.getReadBases(), maxAlignDistance);
    }

    private static int getMaxRealignDistance(final VariantReadContext readContext, final SAMRecord record)
    {
        // returns the larger of (read indel length total + largest core) or (max repeat size  = 5), + 1
        int indelLength = record.getCigar().getCigarElements().stream()
                .filter(x -> x.getOperator() == I || x.getOperator() == D).mapToInt(x -> x.getLength()).sum();

        return max(indelLength + max(readContext.leftCoreLength(), readContext.rightCoreLength()), Realignment.MAX_REPEAT_SIZE) + 1;
    }



    // REALIGN: old logic to be removed
    public static final int MAX_REPEAT_SIZE = 5;
    private static final Repeat NO_REPEAT = new Repeat(0, 0);
    private static final int MIN_REPEAT_COUNT = 4;

    /*
    public static RealignedContext realignedAroundIndex(
            final ReadContext readContext, final int otherBaseIndex, final byte[] otherBases, int maxSize)
    {
        int baseStartIndex = readContext.readBasesLeftFlankIndex();
        int baseEndIndex = readContext.readBasesRightFlankIndex();

        int leftOffset = readContext.readBasesPositionIndex() - baseStartIndex;
        int otherStartIndex = otherBaseIndex - leftOffset;

        return realigned(baseStartIndex, baseEndIndex, readContext.readBases(), otherStartIndex, otherBases, maxSize);
    }
    */

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

        int matchingBasesShortened = matchingBasesFromLeft(baseNextIndex + repeatLength, baseEndIndex, bases, otherNextIndex, otherBases);
        if(matchingBasesShortened > 0 && matchingBases + matchingBasesShortened == exactLength - repeatLength)
        {
            return new RealignedContext(
                    SHORTENED, matchingBasesShortened, otherNextIndex, repeat.RepeatCount, repeatLength, otherIndex, matchingBases);
        }

        int matchingBasesLengthened = matchingBasesFromLeft(baseNextIndex - repeatLength, baseEndIndex, bases, otherNextIndex, otherBases);
        if(matchingBasesLengthened > 0 && matchingBases + matchingBasesLengthened == exactLength + repeatLength)
        {
            return new RealignedContext(
                    LENGTHENED, matchingBasesLengthened, otherNextIndex, repeat.RepeatCount + 1, repeatLength, otherIndex, matchingBases);
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
            int repeats = backwardRepeats(index - i, i, bases) + 1;

            if(repeats >= MIN_REPEAT_COUNT)
                return new Repeat(i, repeats);
        }

        return NO_REPEAT;
    }

    private static int backwardRepeats(int index, int repeatLength, final byte[] readSequence)
    {
        for(int count = 1; ; count++)
        {
            if(!match(index, repeatLength, index - count * repeatLength, readSequence))
            {
                return count - 1;
            }
        }
    }

    public static boolean match(int repeatIndex, int repeatLength, int readIndex, byte[] readSequence)
    {
        return matchingBasesFromLeft(repeatIndex, repeatLength, readIndex, readSequence) == repeatLength;
    }

    private static int matchingBasesFromLeft(int repeatIndex, int repeatLength, int readIndex, byte[] readSequence)
    {
        for(int i = 0; i < repeatLength; i++)
        {
            if(outOfBounds(repeatIndex + i, readSequence) || outOfBounds(readIndex + i, readSequence)
                    || readSequence[repeatIndex + i] != readSequence[readIndex + i])
            {
                return i;
            }
        }

        return repeatLength;
    }

    private static boolean outOfBounds(int index, byte[] sequence)
    {
        return index < 0 || index >= sequence.length;
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

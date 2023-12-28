package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.esvee.read.ReadUtils.getAvgBaseQuality;
import static com.hartwig.hmftools.esvee.read.ReadUtils.isDiscordant;
import static com.hartwig.hmftools.esvee.util.CommonUtils.reverseBytes;

import java.util.Arrays;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.common.Direction;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;

public class ReadRescue
{
    private final RefGenomeInterface mRef;

    public ReadRescue(final RefGenomeInterface ref) {
        mRef = ref;
    }

    public Read rescueRead(final Read read)
    {
        if(read.isUnmapped() || isDiscordant(read) || read.getMappingQuality() < 60)
            return read; // We attempt repair based on the reference genome. If we're not well mapped, this is a terrible idea.

        final Direction direction = read.positiveStrand() ? Direction.FORWARDS : Direction.REVERSE;

        byte repeatBase = 'X';
        int repeatCount = 0;
        for(int i = 0; i < read.getBases().length; i++)
        {
            final int index = direction == Direction.FORWARDS ? i : read.getBases().length - i - 1;
            final byte base = read.getBases()[index];
            if(base == repeatBase)
            {
                repeatCount++;
                continue;
            }

            if(repeatCount >= 10)
            {
                final int averageQuality = direction == Direction.FORWARDS
                        ? getAvgBaseQuality(read, i + 1, read.getBases().length - i + 1)
                        : getAvgBaseQuality(read, 1, index + 1);
                if(averageQuality > 25)
                    continue;

                final Read repaired = tryRescueRead(read, repeatBase, index, direction);
                if(repaired != read)
                    return repaired;
            }
            repeatCount = 0;
            repeatBase = base;
        }

        return read;
    }

    private Read tryRescueRead(final Read read, final byte repeatBase, final int attemptIndex, final Direction direction)
    {
        @Nullable final byte[] newQuals = direction == Direction.FORWARDS
                ? tryRescueReadForwards(read, repeatBase, attemptIndex)
                : tryRescueReadBackwards(read, repeatBase, attemptIndex);
        if(newQuals == null)
            return read;

        /*
        final T clone = (T) read.copyRecord();
        clone.setBases(read.getBases(), newQuals);
        return clone;
         */

        return read;
    }

    @Nullable
    private byte[] tryRescueReadForwards(final Read read, final byte repeatBase, final int attemptIndex)
    {
        final int referencePosition = referencePositionFromRecordIndex(read, attemptIndex + 1);
        final int referenceStartPosition = referencePosition - 5;
        final int referenceEndPosition = referencePosition + 50 + read.getLength();
        if(referenceStartPosition <= 1 || referenceEndPosition >= mRef.getChromosomeLength(read.getChromosome()))
            return null;

        final byte[] referenceBases = mRef.getBases(read.getChromosome(), referenceStartPosition, referenceEndPosition);

        return tryRescueRead(read.getBases(), read.getBaseQuality(), attemptIndex, referenceBases, repeatBase);
    }

    @Nullable
    private byte[] tryRescueReadBackwards(final Read read, final byte repeatBase, final int attemptIndex)
    {
        final int referencePosition = referencePositionFromRecordIndex(read, attemptIndex + 1);
        final int referenceStartPosition = referencePosition - 50 - read.getLength();
        final int referenceEndPosition = referencePosition + 5;
        if(referenceStartPosition <= 1 || referenceEndPosition >= mRef.getChromosomeLength(read.getChromosome()))
            return null;

        final byte[] referenceBases = mRef.getBases(read.getChromosome(), referenceStartPosition, referenceEndPosition);
        final byte[] reversedReferenceBases = reverseBytes(referenceBases);

        final byte[] reversedBases = reverseBytes(read.getBases());
        final byte[] reversedQuals = reverseBytes(read.getBaseQuality());
        final int reversedAttemptIndex = read.getLength() - attemptIndex;
        return reverseBytes(tryRescueRead(reversedBases, reversedQuals, reversedAttemptIndex, reversedReferenceBases, repeatBase));
    }

    private static int referencePositionFromRecordIndex(final Read read, final int desiredReadPosition)
    {
        int readPosition = 1;
        int referencePosition = read.getAlignmentStart();

        for(CigarElement element : read.getCigar().getCigarElements())
        {
            final int nextRefPosition =
                    element.getOperator().consumesReferenceBases() ? referencePosition + element.getLength() : referencePosition;
            final int nextReadPosition = element.getOperator().consumesReadBases() ? readPosition + element.getLength() : readPosition;

            if(nextReadPosition == desiredReadPosition)
                return nextRefPosition;
            else if(nextReadPosition > desiredReadPosition)
            {
                if(element.getOperator().consumesReferenceBases())
                {
                    final int delta = desiredReadPosition - readPosition;
                    return referencePosition + delta;
                }
                else
                    return referencePosition;
            }

            readPosition = nextReadPosition;
            referencePosition = nextRefPosition;
        }

        throw new IllegalStateException("Passed end of read");
    }

    @Nullable
    private byte[] tryRescueRead(
            final byte[] recordBases, final byte[] recordQuals, final int attemptIndex, final byte[] referenceBases, final byte repeatBase)
    {
        // Try to line up recordBases with referenceBases by skipping reference bases
        int bestAgreeCount = 0;
        int bestAgreeReferenceSkip = 0;
        int bestAgreeBaseIndex = 0;
        for(int referenceSkip = 0; referenceSkip < 50; referenceSkip++)
        {
            int agreeCount = 0;
            int strongAgreeNonRepeatCount = 0;
            int strongAgreeNonCCount = 0;
            int strongAgreeNonGCount = 0;
            int lastDisagreeIndex = attemptIndex;
            for(int baseIndex = attemptIndex; baseIndex < recordBases.length; baseIndex++)
            {
                final byte recordBase = recordBases[baseIndex];
                final byte referenceBase = referenceBases[baseIndex - attemptIndex + referenceSkip];
                if(recordBase == referenceBase)
                {
                    agreeCount++;
                    if(recordQuals[baseIndex] >= 25)
                    {
                        if(recordBase != repeatBase)
                            strongAgreeNonRepeatCount++;
                        if(recordBase != 'C')
                            strongAgreeNonCCount++;
                        if(recordBase != 'G')
                            strongAgreeNonGCount++;
                    }
                }
                else if(recordQuals[baseIndex] >= 35)
                {
                    agreeCount = 0;
                    strongAgreeNonCCount = strongAgreeNonGCount = 0;
                    strongAgreeNonRepeatCount = 0;
                    lastDisagreeIndex = baseIndex;
                }
            }

            if(strongAgreeNonRepeatCount >= 4 && strongAgreeNonCCount >= 4 && strongAgreeNonGCount >= 4 && agreeCount > bestAgreeCount)
            {
                bestAgreeCount = agreeCount;
                bestAgreeReferenceSkip = referenceSkip;
                bestAgreeBaseIndex = lastDisagreeIndex;
            }
        }
        if(bestAgreeCount < 16)
            return null;

        // Boost quals in this record.
        final byte[] newQuals = Arrays.copyOf(recordQuals, recordQuals.length);
        for(int baseIndex = bestAgreeBaseIndex; baseIndex < recordBases.length; baseIndex++)
            if(recordBases[baseIndex] == referenceBases[baseIndex - attemptIndex + bestAgreeReferenceSkip])
                newQuals[baseIndex] = (byte) Math.max(newQuals[baseIndex], 37);

        return newQuals;
    }
}

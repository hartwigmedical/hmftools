package com.hartwig.hmftools.svassembly;

import java.util.Arrays;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.sv.Direction;
import com.hartwig.hmftools.svassembly.models.IRecord;
import com.hartwig.hmftools.svassembly.models.MutableRecord;
import com.hartwig.hmftools.svassembly.util.ArrayUtils;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;

public class ReadRescue
{
    private final RefGenomeInterface mRef;

    public ReadRescue(final RefGenomeInterface ref) {
        mRef = ref;
    }

    public <T extends MutableRecord> T rescueRead(final T read)
    {
        if (read.isUnmapped() || read.isDiscordant(1000) || read.getMappingQuality() < 60)
            return read; // We attempt repair based on the reference genome. If we're not well mapped, this is a terrible idea.

        final Direction direction = read.isPositiveStrand() ? Direction.FORWARDS : Direction.REVERSE;

        byte repeatBase = 'X';
        int repeatCount = 0;
        for (int i = 0; i < read.getBases().length; i++)
        {
            final int index = direction == Direction.FORWARDS ? i : read.getBases().length - i - 1;
            final byte base = read.getBases()[index];
            if (base == repeatBase)
            {
                repeatCount++;
                continue;
            }

            if (repeatCount >= 10)
            {
                final int averageQuality = direction == Direction.FORWARDS
                        ? read.getAvgBaseQuality(i + 1, read.getBases().length - i + 1)
                        : read.getAvgBaseQuality(1, index + 1);
                if (averageQuality > 25)
                    continue;

                final T repaired = tryRescueRead(read, repeatBase, index, direction);
                if (repaired != read)
                    return repaired;
            }
            repeatCount = 0;
            repeatBase = base;
        }

        return read;
    }

    private <T extends MutableRecord> T tryRescueRead(final T record, final byte repeatBase, final int attemptIndex, final Direction direction)
    {
        @Nullable final byte[] newQuals = direction == Direction.FORWARDS
                ? tryRescueReadForwards(record, repeatBase, attemptIndex)
                : tryRescueReadBackwards(record, repeatBase, attemptIndex);
        if (newQuals == null)
            return record;

        //noinspection unchecked
        final T clone = (T) record.copy();
        clone.setBases(record.getBases(), newQuals);
        return clone;
    }

    @Nullable
    private byte[] tryRescueReadForwards(final IRecord record, final byte repeatBase, final int attemptIndex)
    {
        final int referencePosition = referencePositionFromRecordIndex(record, attemptIndex + 1);
        final int referenceStartPosition = referencePosition - 5;
        final int referenceEndPosition = referencePosition + 50 + record.getLength();
        if (referenceStartPosition <= 1 || referenceEndPosition >= mRef.getChromosomeLength(record.getChromosome()))
            return null;

        final byte[] referenceBases = mRef.getBases(record.getChromosome(), referenceStartPosition, referenceEndPosition);

        return tryRescueRead(record.getBases(), record.getBaseQuality(), attemptIndex, referenceBases, repeatBase);
    }

    @Nullable
    private byte[] tryRescueReadBackwards(final IRecord record, final byte repeatBase, final int attemptIndex)
    {
        final int referencePosition = referencePositionFromRecordIndex(record, attemptIndex + 1);
        final int referenceStartPosition = referencePosition - 50 - record.getLength();
        final int referenceEndPosition = referencePosition + 5;
        if (referenceStartPosition <= 1 || referenceEndPosition >= mRef.getChromosomeLength(record.getChromosome()))
            return null;

        final byte[] referenceBases = mRef.getBases(record.getChromosome(), referenceStartPosition, referenceEndPosition);
        final byte[] reversedReferenceBases = ArrayUtils.reverse(referenceBases);

        final byte[] reversedBases = ArrayUtils.reverse(record.getBases());
        final byte[] reversedQuals = ArrayUtils.reverse(record.getBaseQuality());
        final int reversedAttemptIndex = record.getLength() - attemptIndex;
        return ArrayUtils.reverse(tryRescueRead(reversedBases, reversedQuals, reversedAttemptIndex, reversedReferenceBases, repeatBase));
    }

    private static int referencePositionFromRecordIndex(final IRecord record, final int desiredReadPosition)
    {
        int readPosition = 1;
        int referencePosition = record.getAlignmentStart();

        for(final CigarElement element : record.getCigar().getCigarElements())
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
    private byte[] tryRescueRead(final byte[] recordBases, final byte[] recordQuals, final int attemptIndex,
            final byte[] referenceBases, final byte repeatBase)
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
                    if (recordQuals[baseIndex] >= 25)
                    {
                        if (recordBase != repeatBase)
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
        for (int baseIndex = bestAgreeBaseIndex; baseIndex < recordBases.length; baseIndex++)
            if(recordBases[baseIndex] == referenceBases[baseIndex - attemptIndex + bestAgreeReferenceSkip])
                newQuals[baseIndex] = (byte) Math.max(newQuals[baseIndex], 37);

        return newQuals;
    }
}

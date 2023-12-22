package com.hartwig.hmftools.svassembly.models;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;

public interface IRecord extends AlignedSequence, Cloneable
{
    boolean isUnmapped();

    String getChromosome();

    int getAlignmentStart();

    int getAlignmentEnd();

    int getUnclippedStart();

    int getUnclippedEnd();

    int getMappingQuality();

    Cigar getCigar();

    boolean isPairedRead();

    boolean isFirstOfPair();

    default boolean isSecondOfPair()
    {
        return isPairedRead() && !isFirstOfPair();
    }

    boolean isGermline();

    String sampleName();

    @Nullable
    String getMateChromosome();

    /** 0 for no alignment */
    int getMateAlignmentStart();

    boolean isMateMapped();

    boolean isMateUnmapped();

    boolean isMatePositiveStrand();

    default boolean isMateOnTheLeft()
    {
        return !isPositiveStrand();
    }

    int impliedFragmentLength();

    boolean isDiscordant(final int discordantPairFragmentLength);

    boolean isPositiveStrand();

    <T> T getAttribute(final String name);

    default int getAvgBaseQuality()
    {
        return getAvgBaseQuality(1, getLength());
    }

    default int getAvgBaseQuality(final int startPosition, final int length)
    {
        final byte[] baseQualities = getBaseQuality();
        final int startIndex = startPosition - 1;
        final int endIndex = Math.min(startIndex + length, baseQualities.length);

        int qualitySum = 0;
        for(int i = startIndex; i < endIndex; i++)
            qualitySum += baseQualities[i];
        return qualitySum / length;
    }

    default int getReadPositionAtReferencePosition(final int pos)
    {
        if(pos <= 0)
            return 0;

        for(final Alignment alignmentBlock : getAlignmentBlocks())
        {
            final int end = alignmentBlock.ReferenceStartPosition + alignmentBlock.Length - 1;
            if(end >= pos)
            {
                if(pos < alignmentBlock.ReferenceStartPosition)
                    return 0;
                else
                    return pos - alignmentBlock.ReferenceStartPosition + alignmentBlock.SequenceStartPosition;
            }
        }
        return 0;
    }

    IRecord copy();
}


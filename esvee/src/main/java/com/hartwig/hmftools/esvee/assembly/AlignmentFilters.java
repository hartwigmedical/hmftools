package com.hartwig.hmftools.esvee.assembly;

import java.util.Arrays;
import java.util.Objects;

import com.hartwig.hmftools.common.sv.Direction;
import com.hartwig.hmftools.esvee.Junction;
import com.hartwig.hmftools.esvee.models.Record;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public enum AlignmentFilters
{
    ;

    public static boolean alignmentCrossesJunction(final Record record, final Junction junction)
    {
        return record.getUnclippedStart() <= junction.position() && junction.position() <= record.getUnclippedEnd();
    }

    public static boolean isRecordAverageQualityAbove(final Record record, final int averageBaseQThreshold)
    {
        return getAvgBaseQuality(record, 1, record.getLength()) >= averageBaseQThreshold;
    }

    private static int getAvgBaseQuality(final Record record, final int readPosition, final int length)
    {
        assert (readPosition > 0);

        final byte[] baseQualities = record.getBaseQuality();
        final int startIndex = readPosition - 1;
        final int endIndex = Math.min(startIndex + length, baseQualities.length);

        int qualitySum = 0;
        for(int i = startIndex; i < endIndex; i++)
            qualitySum += baseQualities[i];
        return qualitySum / length;
    }

    public static boolean isRecordAverageQualityPastJunctionAbove(final Record record, final Junction junction, final int averageBaseQThreshold)
    {
        final int startIndex;
        final int endIndex;
        if(junction.orientation() == Direction.FORWARDS)
        {
            startIndex = junction.position() - record.getUnclippedStart();
            endIndex = record.getLength();
        }
        else
        {
            startIndex = 0;
            endIndex = record.getLength() - (record.getUnclippedEnd() - junction.position());
        }
        if(startIndex == endIndex)
            return true;

        final int averageQuality = getAvgBaseQuality(record, startIndex + 1, endIndex - startIndex);
        return averageQuality > averageBaseQThreshold;
    }

    private static final int ALIGNMENT_SOFT_CLIP_JUNCTION_TOLERANCE = 2;

    public static boolean recordSoftClipsNearJunction(final Record record, final Junction junction)
    {
        // Must start/end with a soft-clip that extends to the junction position
        final CigarElement element = junction.orientation() == Direction.FORWARDS
                ? record.getCigar().getLastCigarElement()
                : record.getCigar().getFirstCigarElement();
        if(element.getOperator() != CigarOperator.S)
            return false;

        final int junctionOffsetFromEnd = junction.orientation() == Direction.FORWARDS
                ? record.getUnclippedEnd() - junction.position()
                : junction.position() - record.getUnclippedStart();
        final int softClipLength = element.getLength();
        return softClipLength + ALIGNMENT_SOFT_CLIP_JUNCTION_TOLERANCE >= junctionOffsetFromEnd;
    }

    public static boolean hasAcceptableMapQ(final Record record, final int threshold)
    {
        return record.getMappingQuality() >= threshold;
    }

    public static boolean isNotBadlyMapped(final Record record)
    {
        return !isBadlyMapped(record);
    }

    private static boolean isBadlyMapped(final Record record)
    {
        if (!record.isDiscordant(1000))
            return false;

        final int mappedSize = record.getCigar().getCigarElements().stream()
                .filter(element -> element.getOperator() == CigarOperator.M)
                .mapToInt(CigarElement::getLength)
                .sum();
        final int indels = record.getCigar().getCigarElements().stream()
                .filter(element -> element.getOperator() == CigarOperator.D || element.getOperator() == CigarOperator.I)
                .mapToInt(CigarElement::getLength)
                .sum();
        final int nm = Objects.requireNonNullElse(record.getAttribute("NM"), indels);

        final int mismatchedBases = nm - indels;
        final int matchedBases = mappedSize - mismatchedBases;
        if (matchedBases < 30)
            return true;
        if (indels > 0 && mismatchedBases >= 3)
            return true;

        final float[] stats = mappedBaseStats(record);
        int countAbove70 = 0;
        int countAbove35 = 0;
        for (final float frequency : stats)
        {
            if(frequency > 0.7f)
                countAbove70++;
            if (frequency > 0.35f)
                countAbove35++;
        }
        //noinspection RedundantIfStatement
        if (countAbove70 >= 1 || countAbove35 >= 2)
            return true;

        return false;
    }

    private static float[] mappedBaseStats(final Record alignment)
    {
        final int[] baseCount = new int[5];

        int readPosition = 1;
        for(final CigarElement element : alignment.getCigar().getCigarElements())
        {
            if (element.getOperator() != CigarOperator.M)
            {
                if (element.getOperator().consumesReadBases())
                    readPosition += element.getLength();
                continue;
            }

            for (int i = 0; i < element.getLength(); i++)
            {
                final byte base = alignment.getBases()[readPosition + i - 1];
                baseCount[baseToIndex(base)]++;
            }
        }

        final float totalBases = Arrays.stream(baseCount).sum();
        final float[] baseFrequency = new float[baseCount.length];
        for (int i = 0; i < baseCount.length; i++)
            baseFrequency[i] = baseCount[i] / totalBases;

        return baseFrequency;
    }

    private static int baseToIndex(final byte base)
    {
        switch(base)
        {
            case 'A':
                return 0;
            case 'T':
                return 1;
            case 'C':
                return 2;
            case 'G':
                return 3;
            default:
                return 4;
        }
    }
}

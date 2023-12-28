package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.esvee.read.ReadUtils.isDiscordant;

import java.util.Arrays;
import java.util.Objects;

import com.hartwig.hmftools.esvee.common.Direction;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.read.Read;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public final class AlignmentFilters
{
    public static boolean alignmentCrossesJunction(final Read read, final Junction junction)
    {
        return read.getUnclippedStart() <= junction.Position && junction.Position <= read.getUnclippedEnd();
    }

    public static boolean isRecordAverageQualityAbove(final Read read, final int averageBaseQThreshold)
    {
        return getAvgBaseQuality(read, 1, read.getLength()) >= averageBaseQThreshold;
    }

    private static int getAvgBaseQuality(final Read read, final int readPosition, final int length)
    {
        assert (readPosition > 0);

        final byte[] baseQualities = read.getBaseQuality();
        final int startIndex = readPosition - 1;
        final int endIndex = Math.min(startIndex + length, baseQualities.length);

        int qualitySum = 0;
        for(int i = startIndex; i < endIndex; i++)
            qualitySum += baseQualities[i];
        return qualitySum / length;
    }

    public static boolean isRecordAverageQualityPastJunctionAbove(final Read read, final Junction junction, final int averageBaseQThreshold)
    {
        final int startIndex;
        final int endIndex;
        if(junction.Orientation == POS_STRAND)
        {
            startIndex = junction.Position - read.getUnclippedStart();
            endIndex = read.getLength();
        }
        else
        {
            startIndex = 0;
            endIndex = read.getLength() - (read.getUnclippedEnd() - junction.position());
        }
        if(startIndex == endIndex)
            return true;

        final int averageQuality = getAvgBaseQuality(read, startIndex + 1, endIndex - startIndex);
        return averageQuality > averageBaseQThreshold;
    }

    private static final int ALIGNMENT_SOFT_CLIP_JUNCTION_TOLERANCE = 2;

    public static boolean recordSoftClipsNearJunction(final Read read, final Junction junction)
    {
        // Must start/end with a soft-clip that extends to the junction position
        final CigarElement element = junction.orientation() == Direction.FORWARDS
                ? read.getCigar().getLastCigarElement()
                : read.getCigar().getFirstCigarElement();
        if(element.getOperator() != CigarOperator.S)
            return false;

        final int junctionOffsetFromEnd = junction.orientation() == Direction.FORWARDS
                ? read.getUnclippedEnd() - junction.position()
                : junction.position() - read.getUnclippedStart();
        final int softClipLength = element.getLength();
        return softClipLength + ALIGNMENT_SOFT_CLIP_JUNCTION_TOLERANCE >= junctionOffsetFromEnd;
    }

    public static boolean hasAcceptableMapQ(final Read read, final int threshold)
    {
        return read.getMappingQuality() >= threshold;
    }

    public static boolean isNotBadlyMapped(final Read read)
    {
        return !isBadlyMapped(read);
    }

    private static boolean isBadlyMapped(final Read read)
    {
        if(!isDiscordant(read))
            return false;

        final int mappedSize = read.getCigar().getCigarElements().stream()
                .filter(element -> element.getOperator() == CigarOperator.M)
                .mapToInt(CigarElement::getLength)
                .sum();

        int indelCount = read.getCigar().getCigarElements().stream()
                .filter(element -> element.getOperator() == CigarOperator.D || element.getOperator() == CigarOperator.I)
                .mapToInt(CigarElement::getLength)
                .sum();

        // CHECK
        Object nmAttribute = read.getAttribute(NUM_MUTATONS_ATTRIBUTE);
        int nmCount = nmAttribute != null ? (int)nmAttribute : indelCount;

        final int mismatchedBases = nmCount - indelCount;
        final int matchedBases = mappedSize - mismatchedBases;

        if(matchedBases < 30)
            return true;

        if(indelCount > 0 && mismatchedBases >= 3)
            return true;

        final float[] stats = mappedBaseStats(read);
        int countAbove70 = 0;
        int countAbove35 = 0;
        for(float frequency : stats)
        {
            if(frequency > 0.7f)
                countAbove70++;
            if(frequency > 0.35f)
                countAbove35++;
        }

        if(countAbove70 >= 1 || countAbove35 >= 2)
            return true;

        return false;
    }

    private static float[] mappedBaseStats(final Read alignment)
    {
        final int[] baseCount = new int[5];

        int readPosition = 1;
        for(CigarElement element : alignment.getCigar().getCigarElements())
        {
            if(element.getOperator() != CigarOperator.M)
            {
                if(element.getOperator().consumesReadBases())
                    readPosition += element.getLength();
                continue;
            }

            for(int i = 0; i < element.getLength(); i++)
            {
                final byte base = alignment.getBases()[readPosition + i - 1];
                baseCount[baseToIndex(base)]++;
            }
        }

        final float totalBases = Arrays.stream(baseCount).sum();
        final float[] baseFrequency = new float[baseCount.length];
        for(int i = 0; i < baseCount.length; i++)
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

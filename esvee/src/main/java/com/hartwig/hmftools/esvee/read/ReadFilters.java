package com.hartwig.hmftools.esvee.read;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.esvee.SvConstants.READ_FILTER_MIN_ALIGNED_BASES;
import static com.hartwig.hmftools.esvee.SvConstants.READ_SOFT_CLIP_JUNCTION_BUFFER;
import static com.hartwig.hmftools.esvee.read.ReadUtils.isDiscordant;

import java.util.Arrays;

import com.hartwig.hmftools.esvee.common.Direction;
import com.hartwig.hmftools.esvee.common.Junction;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public final class ReadFilters
{
    public static boolean alignmentCrossesJunction(final Read read, final Junction junction)
    {
        return positionWithin(junction.Position, read.unclippedStart(), read.unclippedEnd());
    }

    public static boolean isRecordAverageQualityAbove(final byte[] baseQualities, final int averageBaseQThreshold)
    {
        int qualitySum = 0;
        for(int i = 0; i < baseQualities.length; i++)
        {
            qualitySum += baseQualities[i];
        }

        double avgBaseQual = qualitySum / (double)baseQualities.length;
        return avgBaseQual >= averageBaseQThreshold;
    }

    private static int getAvgBaseQuality(final Read read, final int readPosition, final int length)
    {
        int startIndex = readPosition - 1;
        int endIndex = min(startIndex + length, read.getBaseQuality().length);

        int qualitySum = 0;
        for(int i = startIndex; i < endIndex; i++)
        {
            qualitySum += read.getBaseQuality()[i];
        }

        return qualitySum / length;
    }

    public static boolean isRecordAverageQualityPastJunctionAbove(final Read read, final Junction junction, final int averageBaseQThreshold)
    {
        int startIndex;
        int endIndex;

        if(junction.Orientation == POS_STRAND)
        {
            startIndex = junction.Position - read.unclippedStart();
            endIndex = read.basesLength();
        }
        else
        {
            startIndex = 0;
            endIndex = read.basesLength() - (read.unclippedEnd() - junction.position());
        }
        if(startIndex == endIndex)
            return true;

        final int averageQuality = getAvgBaseQuality(read, startIndex + 1, endIndex - startIndex);
        return averageQuality > averageBaseQThreshold;
    }

    public static boolean recordSoftClipsNearJunction(final Read read, final Junction junction)
    {
        // previous method was almost certainly invalid
        if(junction.isForward())
        {
            if(!read.isRightClipped())
                return false;

            int junctionAlignmentDiff = abs(read.alignmentEnd() - junction.Position);
            return junctionAlignmentDiff <= READ_SOFT_CLIP_JUNCTION_BUFFER && read.unclippedEnd() > junction.Position;
        }
        else
        {
            if(!read.isLeftClipped())
                return false;

            int junctionAlignmentDiff = abs(read.alignmentStart() - junction.Position);
            return junctionAlignmentDiff <= READ_SOFT_CLIP_JUNCTION_BUFFER && read.unclippedStart() < junction.Position;
        }
    }

    public static boolean hasAcceptableMapQ(final Read read, final int threshold)
    {
        return read.mappingQuality() >= threshold;
    }

    public static boolean isNotBadlyMapped(final Read read)
    {
        return !isBadlyMapped(read);
    }

    public static boolean isBadlyMapped(final Read read)
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

        if(matchedBases < READ_FILTER_MIN_ALIGNED_BASES)
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

package com.hartwig.hmftools.esvee.old;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.esvee.SvConstants.BAM_READ_JUNCTION_BUFFER;
import static com.hartwig.hmftools.esvee.SvConstants.READ_SOFT_CLIP_JUNCTION_BUFFER;

import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.read.Read;

public final class ReadFilters
{
    public static boolean isCandidateJunctionRead(final Read read, final Junction junction)
    {
        int junctionBoundaryStart = junction.isForward() ? junction.Position - BAM_READ_JUNCTION_BUFFER : junction.Position;
        int junctionBoundaryEnd = junction.isForward() ? junction.Position : junction.Position + BAM_READ_JUNCTION_BUFFER;

        return positionsOverlap(junctionBoundaryStart, junctionBoundaryEnd, read.unclippedStart(), read.unclippedEnd());
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
}

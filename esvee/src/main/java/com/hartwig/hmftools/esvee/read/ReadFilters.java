package com.hartwig.hmftools.esvee.read;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.SvConstants.AVG_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.SvConstants.READ_SOFT_CLIP_JUNCTION_BUFFER;

import com.hartwig.hmftools.esvee.common.Junction;

public final class ReadFilters
{
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

    public static boolean isAboveBaseQualAvgThreshold(final byte[] baseQualities)
    {
        int qualitySum = 0;
        for(int i = 0; i < baseQualities.length; i++)
        {
            qualitySum += baseQualities[i];
        }

        double avgBaseQual = qualitySum / (double)baseQualities.length;
        return avgBaseQual >= AVG_BASE_QUAL_THRESHOLD;
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
}

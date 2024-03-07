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

    public static boolean recordSoftClipsAndCrossesJunction(final Read read, final Junction junction)
    {
        if(junction.isForward())
        {
            return read.isRightClipped() && read.unclippedEnd() > junction.Position;

            //int junctionAlignmentDiff = abs(read.alignmentEnd() - junction.Position);
            //return junctionAlignmentDiff <= READ_SOFT_CLIP_JUNCTION_BUFFER && read.unclippedEnd() > junction.Position;
        }
        else
        {
            return read.isLeftClipped() && read.unclippedStart() < junction.Position;
        }
    }

    public static boolean recordSoftClipsAtJunction(final Read read, final Junction junction)
    {
        if(junction.isForward())
        {
            return read.isRightClipped() && read.alignmentEnd() == junction.Position;
        }
        else
        {
            return read.isLeftClipped() && read.alignmentStart() == junction.Position;
        }
    }

}

package com.hartwig.hmftools.esvee.read;

import static java.lang.Math.max;
import static java.lang.Math.min;

import com.hartwig.hmftools.esvee.types.Junction;

public final class ReadFilters
{
    public static boolean recordSoftClipsAndCrossesJunction(final Read read, final Junction junction)
    {
        if(junction.isForward())
        {
            return read.isRightClipped() && read.unclippedEnd() > junction.Position;
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

    public static int readJunctionExtensionLength(final Read read, final Junction junction)
    {
        if(junction.isForward())
        {
            return read.isRightClipped() ? max(read.unclippedEnd() - junction.Position, 0) : 0;
        }
        else
        {
            return read.isLeftClipped() ? max(junction.Position - read.unclippedStart(), 0) : 0;
        }
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

    public static boolean isAboveBaseQualAvgThreshold(final byte[] baseQualities, final int threshold)
    {
        int qualitySum = 0;
        for(int i = 0; i < baseQualities.length; i++)
        {
            qualitySum += baseQualities[i];
        }

        double avgBaseQual = qualitySum / (double)baseQualities.length;
        return avgBaseQual >= threshold;
    }
}

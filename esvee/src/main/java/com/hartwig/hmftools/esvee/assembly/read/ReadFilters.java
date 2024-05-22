package com.hartwig.hmftools.esvee.assembly.read;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MAX_NON_SOFT_CLIP_OVERLAP;

import com.hartwig.hmftools.esvee.assembly.types.Junction;

public final class ReadFilters
{
    public static boolean recordSoftClipsAndCrossesJunction(final Read read, final Junction junction)
    {
        int junctionOverlap = 0;

        if(junction.isForward())
        {
            if(read.isRightClipped())
                return read.maxUnclippedEnd() > junction.Position;

            // consider alignments with small junction overlap, which ought to have mismatched bases
            junctionOverlap = read.alignmentEnd() - junction.Position;
        }
        else
        {
            if(read.isLeftClipped())
                return read.minUnclippedStart() < junction.Position;

            junctionOverlap = junction.Position - read.alignmentStart();
        }

        return read.snvCount() > 0 && junctionOverlap > 0 && junctionOverlap <= PRIMARY_ASSEMBLY_MAX_NON_SOFT_CLIP_OVERLAP;
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
            return read.isRightClipped() ? max(read.maxUnclippedEnd() - junction.Position, 0) : 0;
        }
        else
        {
            return read.isLeftClipped() ? max(junction.Position - read.minUnclippedStart(), 0) : 0;
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

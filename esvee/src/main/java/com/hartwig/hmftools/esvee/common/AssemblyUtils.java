package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.read.Read;

public final class AssemblyUtils
{
    public static JunctionAssembly buildFromJunctionReads(final Junction junction, final List<Read> reads, boolean checkMismatches)
    {
        JunctionAssembly junctionSequence = buildFromJunction(junction, reads);

        for(Read read : reads)
        {
            if(read == junctionSequence.initialRead())
                continue;

            junctionSequence.addRead(read, checkMismatches);
        }

        return junctionSequence;
    }

    protected static JunctionAssembly buildFromJunction(final Junction junction, final List<Read> reads)
    {
        int minAlignedPosition = junction.Position;
        int maxAlignedPosition = junction.Position;

        Read maxJunctionBaseQualRead = null;
        int maxJunctionBaseQualTotal = 0;

        for(Read read : reads)
        {
            if(junction.isForward())
            {
                maxAlignedPosition = max(maxAlignedPosition, read.unclippedEnd());
            }
            else
            {
                minAlignedPosition = min(minAlignedPosition, read.unclippedStart());
            }

            int junctionBaseQualTotal = readQualFromJunction(read, junction);

            if(junctionBaseQualTotal > maxJunctionBaseQualTotal)
            {
                maxJunctionBaseQualTotal = junctionBaseQualTotal;
                maxJunctionBaseQualRead = read;
            }
        }

        return new JunctionAssembly(junction, maxJunctionBaseQualRead, minAlignedPosition,  maxAlignedPosition);
    }

    protected static int readQualFromJunction(final Read read, final Junction junction)
    {
        int junctionReadIndex = read.getReadIndexAtReferencePosition(junction.Position, true);

        int readIndexStart;
        int readIndexEnd;

        if(junction.direction() == Direction.FORWARDS)
        {
            readIndexStart = junctionReadIndex;
            readIndexEnd = read.basesLength() - 1;
        }
        else
        {
            readIndexStart = 0;
            readIndexEnd = junctionReadIndex;
        }

        int baseQualTotal = 0;

        for(int i = readIndexStart; i <= readIndexEnd; ++i)
        {
            baseQualTotal += read.getBaseQuality()[i];
        }

        return baseQualTotal;
    }

    public static boolean basesMatch(
            final byte first, final byte second, final byte firstQual, final byte secondQual, final int lowQualThreshold)
    {
        return first == second || firstQual < lowQualThreshold || secondQual < lowQualThreshold;
    }
}

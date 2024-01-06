package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.List;

import com.hartwig.hmftools.esvee.read.Read;

public final class AssemblyUtils
{
    public static AssemblySequence buildFromJunctionReads(final Junction junction, final List<Read> reads, boolean checkMismatches)
    {
        AssemblySequence junctionSequence = buildFromJunction(junction, reads);

        for(Read read : reads)
        {
            if(read == junctionSequence.initialRead())
                continue;

            junctionSequence.addRead(read, checkMismatches);
        }

        return junctionSequence;
    }

    protected static AssemblySequence buildFromJunction(final Junction junction, final List<Read> reads)
    {
        int minAlignedPosition = junction.Position;
        int maxAlignedPosition = junction.Position;

        Read maxJunctionBaseQualRead = null;
        int maxJunctionBaseQualTotal = 0;

        for(Read read : reads)
        {
            if(junction.direction() == Direction.FORWARDS)
            {
                maxAlignedPosition = max(maxAlignedPosition, read.getUnclippedEnd());
            }
            else
            {
                minAlignedPosition = min(minAlignedPosition, read.getUnclippedStart());
            }

            int junctionBaseQualTotal = readQualFromJunction(read, junction);

            if(junctionBaseQualTotal > maxJunctionBaseQualTotal)
            {
                maxJunctionBaseQualTotal = junctionBaseQualTotal;
                maxJunctionBaseQualRead = read;
            }
        }

        return new AssemblySequence(junction, maxJunctionBaseQualRead, minAlignedPosition,  maxAlignedPosition);
    }

    protected static int readQualFromJunction(final Read read, final Junction junction)
    {
        int junctionReadIndex = read.getReadIndexAtReferencePosition(junction.Position, true);

        int readIndexStart;
        int readIndexEnd;

        if(junction.direction() == Direction.FORWARDS)
        {
            readIndexStart = junctionReadIndex;
            readIndexEnd = read.getLength() - 1;
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

    public static boolean purgeLowSupport(final AssemblySequence assemblySequence, int minReadCount, int minTotalQual)
    {
        boolean hasSupportedMismatches = false;

        for(BaseMismatches baseMismatches : assemblySequence.mismatches().baseMismatches())
        {
            for(int j = 0; j < baseMismatches.Mismatches.length; ++j)
            {
                if(baseMismatches.Mismatches[j] == null)
                    continue;

                BaseMismatch mismatch = baseMismatches.Mismatches[j];

                if(mismatch.Reads.size() >= minReadCount && mismatch.QualTotal >= minTotalQual)
                {
                    hasSupportedMismatches = true;
                }
                else
                {
                    // remove the mismatch
                    baseMismatches.Mismatches[j].Reads.clear();
                    baseMismatches.Mismatches[j] = null;
                }
            }
        }

        return hasSupportedMismatches;
    }

}

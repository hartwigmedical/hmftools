package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.List;

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
        // find the longest extension out from the junction into the unmapped region ie opposite to the orientation
        int minAlignedPosition = junction.Position;
        int maxAlignedPosition = junction.Position;

        Read maxJunctionBaseQualRead = null;
        int maxJunctionBaseQualTotal = 0;

        int maxDistanceFromJunction = 0;

        for(Read read : reads)
        {
            int readJunctionIndex = read.getReadIndexAtReferencePosition(junction.Position, true);

            // calculate how many bases beyond the junction the read extends
            // for positive orientations, if read length is 10, and junction index is at 6, then extends with indices 7-9 ie 3
            // for negative orientations, if read length is 10, and junction index is at 4, then extends with indices 0-3 ie 4
            int extensionDistance = junction.isForward() ? read.basesLength() - readJunctionIndex - 1 : readJunctionIndex;

            maxDistanceFromJunction = max(maxDistanceFromJunction, extensionDistance);

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

        return new JunctionAssembly(junction, maxJunctionBaseQualRead, maxDistanceFromJunction, minAlignedPosition,  maxAlignedPosition);
    }

    public static void expandReferenceBases(final JunctionAssembly assembly)
    {
        // find the longest length of reference bases extending back from the junction
        int minAlignedPosition = assembly.minAlignedPosition();
        int maxAlignedPosition = assembly.maxAlignedPosition();

        int junctionPosition = assembly.junction().Position;
        boolean junctionIsForward = assembly.junction().isForward();

        int maxDistanceFromJunction = 0;

        AssemblySupport minNmSupport = null;
        int minNmSupportMaxDistance = 0;

        for(AssemblySupport support : assembly.support())
        {
            int readJunctionIndex = support.read().getReadIndexAtReferencePosition(junctionPosition, true);

            // for positive orientations, if read length is 10, and junction index is at 4, then extends with indices 0-3 ie 4
            // for negative orientations, if read length is 10, and junction index is at 6, then extends with indices 7-9 ie 4
            int extensionDistance = junctionIsForward ? readJunctionIndex : support.read().basesLength() - readJunctionIndex - 1;

            maxDistanceFromJunction = max(maxDistanceFromJunction, extensionDistance);

            if(junctionIsForward)
            {
                minAlignedPosition = min(minAlignedPosition, support.read().unclippedStart());
            }
            else
            {
                maxAlignedPosition = max(maxAlignedPosition, support.read().unclippedEnd());
            }

            if(minNmSupport == null
            || support.read().numberOfEvents() < minNmSupport.read().numberOfEvents()
            || (support.read().numberOfEvents() == minNmSupport.read().numberOfEvents() && extensionDistance < minNmSupportMaxDistance))
            {
                minNmSupport = support;
                minNmSupportMaxDistance = extensionDistance;
            }
        }

        assembly.extendBases(maxDistanceFromJunction, minAlignedPosition, maxAlignedPosition);

        // order by NM to favour the ref where possible
        if(minNmSupport != null)
        {
            assembly.extendReadSupport(minNmSupport.read(), minNmSupport);
        }

        for(AssemblySupport support : assembly.support())
        {
            if(support == minNmSupport)
                continue;

            assembly.extendReadSupport(support.read(), support);
        }
    }

    protected static int readQualFromJunction(final Read read, final Junction junction)
    {
        int readJunctionIndex = read.getReadIndexAtReferencePosition(junction.Position, true);

        int readIndexStart;
        int readIndexEnd;

        if(junction.direction() == Direction.FORWARDS)
        {
            readIndexStart = readJunctionIndex;
            readIndexEnd = read.basesLength() - 1;
        }
        else
        {
            readIndexStart = 0;
            readIndexEnd = readJunctionIndex;
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

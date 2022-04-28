package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.isofox.common.CommonUtils.deriveCommonRegions;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.isofox.common.ReadRecord;

public class LocalJunctionData
{
    public final int[] JunctionPositions;
    public final byte[] JunctionOrientations;
    public int MatchCount;
    public final int[] MaxSplitLengths;

    public LocalJunctionData(final int[] junctionPositions, final byte[] junctionOrientations)
    {
        JunctionPositions = junctionPositions;
        JunctionOrientations = junctionOrientations;
        MaxSplitLengths = new int[SE_PAIR];
        MatchCount = 0;
    }

    public String toString()
    {
        return String.format("junc(%d - %d) matched(%d) maxSplits(%d - %d)",
                JunctionPositions[SE_START], JunctionPositions[SE_END], MatchCount,
                MaxSplitLengths[SE_START], MaxSplitLengths[SE_END]);
    }

    public static void setMaxSplitMappedLength(
            int seIndex, final List<ReadRecord> reads, final int[] junctPositions, final byte[] junctOrientations, final int[] maxSplitLengths)
    {
        // find the longest section mapped across the junction
        final List<ReadRecord> matchingReads = reads.stream()
                .filter(x -> positionWithin(junctPositions[seIndex], x.PosStart, x.PosEnd)).collect(Collectors.toList());

        if(matchingReads.isEmpty()) // can occur with the fragments from a fusion merged in due to homology
            return;

        List<int[]> mappedCoords;

        if(matchingReads.size() == 1)
        {
            mappedCoords = matchingReads.get(0).getMappedRegionCoords(false);
        }
        else
        {
            mappedCoords = deriveCommonRegions(
                    matchingReads.get(0).getMappedRegionCoords(false), matchingReads.get(1).getMappedRegionCoords(false));
        }

        int mappedBases = 0;

        for(int[] coord : mappedCoords)
        {
            if(junctOrientations[seIndex] == NEG_ORIENT)
            {
                if(coord[SE_END] < junctPositions[seIndex])
                    continue;

                mappedBases += coord[SE_END] - max(junctPositions[seIndex], coord[SE_START]) + 1;
            }
            else
            {
                if(coord[SE_START] > junctPositions[seIndex])
                    break;

                mappedBases += min(junctPositions[seIndex], coord[SE_END]) - coord[SE_START] + 1;
            }
        }

        maxSplitLengths[seIndex] = max(mappedBases, maxSplitLengths[seIndex]);
    }

}

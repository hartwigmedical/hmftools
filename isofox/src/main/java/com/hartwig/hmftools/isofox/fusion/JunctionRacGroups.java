package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.SOFT_CLIP_JUNC_BUFFER;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.softClippedReadSupportsJunction;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class JunctionRacGroups
{
    // keyed by junction position
    private final Map<Integer,List<ReadGroup>> mPosJunctionGroups;
    private final Map<Integer,List<ReadGroup>> mNegJunctionGroups;

    public JunctionRacGroups()
    {
        mPosJunctionGroups = Maps.newHashMap();
        mNegJunctionGroups = Maps.newHashMap();
    }

    public Map<Integer,List<ReadGroup>> getJunctionGroups(byte orientation)
    {
        return orientation == POS_ORIENT ? mPosJunctionGroups : mNegJunctionGroups;
    }

    public void clear()
    {
        mPosJunctionGroups.clear();
        mNegJunctionGroups.clear();
    }

    public void purgeGroups(int position)
    {
        List<Integer> pastJuncPositions = mPosJunctionGroups.keySet().stream()
                .filter(x -> x < position).collect(Collectors.toList());

        pastJuncPositions.forEach(x -> mPosJunctionGroups.remove(x));

        pastJuncPositions = mNegJunctionGroups.keySet().stream()
                .filter(x -> x < position).collect(Collectors.toList());

        pastJuncPositions.forEach(x -> mNegJunctionGroups.remove(x));
    }

    public void addJunction(int position, byte orientation)
    {
        Map<Integer,List<ReadGroup>> junctionGroups = orientation == POS_ORIENT ? mPosJunctionGroups : mNegJunctionGroups;

        if(!junctionGroups.containsKey(position))
        {
            junctionGroups.put(position, Lists.newArrayList());
        }
    }

    public void checkAddCandidateGroup(final ReadGroup readGroup)
    {
        for(int i = 0; i <= 1; ++i)
        {
            byte juncOrientation = (i == 0) ? POS_ORIENT : NEG_ORIENT;
            Map<Integer,List<ReadGroup>> junctionGroups = (i == 0) ? mPosJunctionGroups : mNegJunctionGroups;
            int seIndex = (i == 0) ? SE_START : SE_END;

            for(Integer juncPosition : junctionGroups.keySet())
            {
                // check that none of the other reads are on the incorrect side of this fusion junction
                if(juncOrientation == POS_ORIENT)
                {
                    if(readGroup.Reads.stream().anyMatch(x -> x.getCoordsBoundary(SE_END) > juncPosition + SOFT_CLIP_JUNC_BUFFER))
                        continue;
                }
                else
                {
                    if(readGroup.Reads.stream().anyMatch(x -> x.getCoordsBoundary(SE_START) < juncPosition - SOFT_CLIP_JUNC_BUFFER))
                        continue;
                }

                if(readGroup.Reads.stream().
                        anyMatch(x -> softClippedReadSupportsJunction(x, seIndex, juncPosition, juncOrientation, null)))
                {
                    junctionGroups.get(juncPosition).add(readGroup);
                }
            }
        }
    }

}

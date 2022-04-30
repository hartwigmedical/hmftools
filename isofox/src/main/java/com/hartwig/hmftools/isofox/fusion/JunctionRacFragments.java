package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.SOFT_CLIP_JUNC_BUFFER;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGNED;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGN_CANDIDATE;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.softClippedReadSupportsJunction;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class JunctionRacFragments
{
    // keyed by junction position
    private final Map<Integer,List<FusionFragment>> mPosJunctionFragments;
    private final Map<Integer,List<FusionFragment>> mNegJunctionFragments;

    private final List<FusionFragment> mRacFragments;

    public JunctionRacFragments()
    {
        mPosJunctionFragments = Maps.newHashMap();
        mNegJunctionFragments = Maps.newHashMap();
        mRacFragments = Lists.newArrayList();
    }

    public Map<Integer,List<FusionFragment>> getJunctionGroups(byte orientation)
    {
        return orientation == POS_ORIENT ? mPosJunctionFragments : mNegJunctionFragments;
    }

    public int fragmentCount() { return mRacFragments.size(); }
    public int assignedFragmentCount() { return (int)mRacFragments.stream().filter(x -> x.type() == REALIGNED).count(); }

    public List<FusionFragment> getRacGroups() { return mRacFragments; }
    public int junctionCount() { return mPosJunctionFragments.size() + mNegJunctionFragments.size(); }

    public List<FusionFragment> getJunctionFragments(byte orientation, int position)
    {
        if(orientation == POS_ORIENT)
            return mPosJunctionFragments.get(position);
        else
            return mNegJunctionFragments.get(position);
    }

    public void clear()
    {
        mPosJunctionFragments.clear();
        mNegJunctionFragments.clear();
        mRacFragments.clear();
    }

    public void addJunction(int position, byte orientation)
    {
        Map<Integer,List<FusionFragment>> junctionGroups = orientation == POS_ORIENT ? mPosJunctionFragments : mNegJunctionFragments;

        if(!junctionGroups.containsKey(position))
        {
            junctionGroups.put(position, Lists.newArrayList());
        }
    }

    public boolean checkAddCandidateGroup(final ReadGroup readGroup)
    {
        boolean fragmentAdded = false;

        for(int i = 0; i <= 1; ++i)
        {
            byte juncOrientation = (i == 0) ? POS_ORIENT : NEG_ORIENT;
            Map<Integer,List<FusionFragment>> junctionGroups = (i == 0) ? mPosJunctionFragments : mNegJunctionFragments;
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
                    FusionFragment fragment = new FusionFragment(readGroup);

                    if(fragment.type() != REALIGN_CANDIDATE)
                        continue;

                    junctionGroups.get(juncPosition).add(fragment);

                    if(!fragmentAdded)
                    {
                        fragmentAdded = true;
                        mRacFragments.add(fragment);
                    }
                }
            }
        }

        return fragmentAdded;
    }

    @Deprecated
    public void purgeGroups(int position)
    {
        List<Integer> pastJuncPositions = mPosJunctionFragments.keySet().stream()
                .filter(x -> x < position).collect(Collectors.toList());

        pastJuncPositions.forEach(x -> mPosJunctionFragments.remove(x));

        pastJuncPositions = mNegJunctionFragments.keySet().stream()
                .filter(x -> x < position).collect(Collectors.toList());

        pastJuncPositions.forEach(x -> mNegJunctionFragments.remove(x));
    }


}

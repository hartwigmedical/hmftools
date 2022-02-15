package com.hartwig.hmftools.sage.phase;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

public class PhasedGroupCollection
{
    private final int mId;
    private final List<PhasedVariantGroup> mGroups;
    private final Map<Integer,List<PhasedVariantGroup>> mGroupsMap; // keyed by min +ve position
    private int mMinPostion;
    private int mMaxPostion;

    public PhasedGroupCollection(int id)
    {
        mId = id;
        mGroups = Lists.newArrayList();
        mGroupsMap = Maps.newHashMap();
        mMinPostion = 0;
        mMaxPostion = 0;
    }

    public int minPosition() { return mMinPostion;}
    public int maxPosition() { return mMaxPostion;}

    public List<PhasedVariantGroup> groups() { return mGroups; }
    public Map<Integer,List<PhasedVariantGroup>> groupsMap() { return mGroupsMap; }

    public void finalise()
    {
        mGroupsMap.values().forEach(x -> mGroups.addAll(x));
        Collections.sort(mGroups, new PhasedVariantGroup.PhasedGroupComparator());
        mGroupsMap.clear();
    }

    public boolean positionsOverlap(int otherMin, int otherMax)
    {
        return BaseRegion.positionsOverlap(mMinPostion, mMaxPostion, otherMin, otherMax);
    }

    public boolean addPhaseVariants(
            int posVarMin, int posVarMax, int nextGroupId,
            final List<ReadContextCounter> posCounters, final List<ReadContextCounter> negCounters)
    {
        List<PhasedVariantGroup> groups = mGroupsMap.get(posVarMin);

        if(groups == null)
        {
            groups = Lists.newArrayList();
            mGroupsMap.put(posVarMin, groups);
        }

        for(PhasedVariantGroup group : groups)
        {
            if(group.exactMatch(posVarMin, posVarMax, posCounters, negCounters))
            {
                group.ReadCount++;
                group.mergeNegatives(negCounters);
                return false;
            }
        }

        PhasedVariantGroup newGroup = new PhasedVariantGroup(nextGroupId, posVarMin, posVarMax, posCounters, negCounters);
        groups.add(newGroup);

        mMinPostion = min(mMinPostion, newGroup.variantMin());
        mMaxPostion = max(mMaxPostion, newGroup.variantMax());
        return true;
    }

    public boolean addPhaseVariantsOld(
            int posVarMin, int posVarMax, int nextGroupId,
            final List<ReadContextCounter> posCounters, final List<ReadContextCounter> negCounters)
    {
        // returns true if a new group was added
        if(mGroups.isEmpty())
        {
            PhasedVariantGroup newGroup = new PhasedVariantGroup(nextGroupId, posVarMin, posVarMax, posCounters, negCounters);
            mGroups.add(newGroup);

            mMinPostion = newGroup.variantMin();
            mMaxPostion = newGroup.variantMax();
            return true;
        }

        int index = 0;

        while(index < mGroups.size())
        {
            PhasedVariantGroup group = mGroups.get(index);

            if(posVarMin < group.posVariantMin())
                break;

            if(group.posVariantMin() == posVarMin)
            {
                // test for an exact match
                if(group.exactMatch(posVarMin, posVarMax, posCounters, negCounters))
                {
                    group.ReadCount++;
                    group.mergeNegatives(negCounters);
                    return false;
                }
            }

            index++;
        }

        if(index > 0 && index < mGroups.size() - 1)
        {
            if(posVarMin < mGroups.get(index - 1).posVariantMin() || posVarMin > mGroups.get(index + 1).posVariantMin())
            {
                SG_LOGGER.error("invalid pos() insertion at index({}) groups({})",
                        posVarMin, index, mGroups.size());
            }
        }

        PhasedVariantGroup newGroup = new PhasedVariantGroup(nextGroupId, posVarMin, posVarMax, posCounters, negCounters);
        mGroups.add(index, newGroup);

        mMinPostion = min(mMinPostion, newGroup.variantMin());
        mMaxPostion = max(mMaxPostion, newGroup.variantMax());
        return true;
    }

    public void merge(final PhasedGroupCollection other)
    {
        mGroupsMap.putAll(other.groupsMap());
        // mGroups.addAll(other.groups());
        // Collections.sort(mGroups, new PhasedVariantGroup.PhasedGroupComparator());
        mMinPostion = min(mMinPostion, other.minPosition());
        mMaxPostion = max(mMaxPostion, other.maxPosition());
    }

    public int groupCount() { return !mGroupsMap.isEmpty() ? mGroupsMap.values().stream().mapToInt(x -> x.size()).sum() : mGroups.size(); }

    public String toString()
    {
        return String.format("id(%d) groups(%d) range(%d - %d)", mId, groupCount(), mMinPostion, mMaxPostion);
    }
}

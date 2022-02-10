package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

public class PhasedGroupCollection
{
    private final List<PhasedVariantGroup> mGroups;
    private int mMinPostion;
    private int mMaxPostion;

    public PhasedGroupCollection()
    {
        mGroups = Lists.newArrayList();
        mMinPostion = 0;
        mMaxPostion = 0;
    }

    public int minPosition() { return mMinPostion;}
    public int maxPosition() { return mMaxPostion;}
    public List<PhasedVariantGroup> groups() { return mGroups; }

    public boolean positionsOverlap(int otherMin, int otherMax)
    {
        return BaseRegion.positionsOverlap(mMinPostion, mMaxPostion, otherMin, otherMax);
    }

    public boolean addPhaseVariants(
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
        mGroups.addAll(other.groups());
        Collections.sort(mGroups, new PhasedVariantGroup.PhasedGroupComparator());
        mMinPostion = min(mMinPostion, other.minPosition());
        mMaxPostion = max(mMaxPostion, other.maxPosition());
    }

    public String toString() { return String.format("groups(%d) range(%d - %d)", mGroups.size(), mMinPostion, mMaxPostion); }
}

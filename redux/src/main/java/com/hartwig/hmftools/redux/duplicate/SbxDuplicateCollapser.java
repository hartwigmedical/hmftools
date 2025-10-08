package com.hartwig.hmftools.redux.duplicate;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.sin;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.collect.MergeUtils.clusterMerger;
import static com.hartwig.hmftools.redux.duplicate.CollapseUtils.DUPLICATE_GROUP_MERGER;
import static com.hartwig.hmftools.redux.duplicate.DuplicateGroup.DUPLICATE_GROUP_COMPARATOR;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.function.BiPredicate;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.redux.common.ReadInfo;

import org.jetbrains.annotations.Nullable;

public class SbxDuplicateCollapser
{
    private final int mMaxDuplicateDistance;

    public SbxDuplicateCollapser(int maxDuplicateDistance)
    {
        mMaxDuplicateDistance = maxDuplicateDistance;
    }

    private class GroupInfo implements Comparable<GroupInfo>
    {
        public final Object Group;
        public final int PositionLower;
        public final int PositionUpper;

        private int mSize;
        private List<GroupInfo> mCollapsedGroups;

        public GroupInfo(final Object group, final int positionLower, final int positionUpper, final int initialSize)
        {
            Group = group;
            PositionLower = positionLower;
            PositionUpper = positionUpper;
            mSize = initialSize;
            mCollapsedGroups = null;
        }

        public int size() { return mSize; }
        public void markCollapsed() { mSize = 0; }
        public boolean collapsed() { return mSize == 0; }

        public boolean hasCollapsedGroups() { return mCollapsedGroups != null; }
        public List<GroupInfo> collapsedGroups() { return mCollapsedGroups; }

        public void addGroup(final GroupInfo group)
        {
            if(mCollapsedGroups == null)
                mCollapsedGroups = Lists.newArrayListWithCapacity(2);

            mCollapsedGroups.add(group);
            mSize += group.size();
        }

        public boolean withinRange(final GroupInfo other)
        {
            return abs(PositionLower - other.PositionLower) + abs(PositionUpper - other.PositionUpper) <= mMaxDuplicateDistance;
        }

        @Override
        public int compareTo(final GroupInfo o)
        {
            int diffLower = PositionLower - o.PositionLower;
            if(diffLower != 0)
            {
                return diffLower;
            }

            return PositionUpper - o.PositionUpper;
        }

        public String toString() { return format("coords(%d - %d) size(%d)", PositionLower, PositionUpper, mSize); }
    }

    public FragmentCoordReads collapseGroups(
            @Nullable final List<DuplicateGroup> duplicateGroups, @Nullable final List<ReadInfo> singleReads)
    {
        Map<String,List<GroupInfo>> keyGroups = Maps.newHashMap();

        int maxCapacity = (singleReads != null ? singleReads.size() : 0) + (duplicateGroups != null ? duplicateGroups.size() : 0);

        if(singleReads != null)
        {
            for(ReadInfo readInfo : singleReads)
            {
                String key = collapseKey(readInfo.coordinates());

                List<GroupInfo> groups = keyGroups.get(key);

                if(groups == null)
                {
                    groups = Lists.newArrayListWithCapacity(maxCapacity);
                    keyGroups.put(key, groups);
                }

                GroupInfo groupInfo = new GroupInfo(
                        readInfo, readInfo.coordinates().PositionLower, readInfo.coordinates().PositionUpper, 1);

                groups.add(groupInfo);
            }
        }

        if(duplicateGroups != null)
        {
            for(DuplicateGroup duplicateGroup : duplicateGroups)
            {
                String key = collapseKey(duplicateGroup.fragmentCoordinates());

                List<GroupInfo> groups = keyGroups.get(key);

                if(groups == null)
                {
                    groups = Lists.newArrayListWithCapacity(maxCapacity);
                    keyGroups.put(key, groups);
                }

                GroupInfo groupInfo = new GroupInfo(
                        duplicateGroup, duplicateGroup.fragmentCoordinates().PositionLower,
                        duplicateGroup.fragmentCoordinates().PositionUpper, duplicateGroup.reads().size());

                groups.add(groupInfo);
            }
        }

        List<DuplicateGroup> finalDuplicateGroups = null;
        List<ReadInfo> finalSingleReads = null;

        for(List<GroupInfo> groups : keyGroups.values())
        {
            collapseGroups(groups);

            for(GroupInfo groupInfo : groups)
            {
                if(groupInfo.collapsed())
                    continue;

                DuplicateGroup duplicateGroup = null;

                if(groupInfo.Group instanceof ReadInfo)
                {
                    ReadInfo readInfo = (ReadInfo)groupInfo.Group;

                    if(!groupInfo.hasCollapsedGroups())
                    {
                        if(finalSingleReads == null)
                            finalSingleReads = Lists.newArrayListWithCapacity(maxCapacity);

                        finalSingleReads.add(readInfo);
                        continue;
                    }

                    duplicateGroup = new DuplicateGroup(null, readInfo.read(), readInfo.coordinates());
                }
                else
                {
                    duplicateGroup = (DuplicateGroup)groupInfo.Group;
                }

                if(groupInfo.hasCollapsedGroups())
                    collapseToDuplicateGroup(duplicateGroup, groupInfo.collapsedGroups());

                if(finalDuplicateGroups == null)
                    finalDuplicateGroups = Lists.newArrayListWithCapacity(maxCapacity);

                finalDuplicateGroups.add(duplicateGroup);
            }
        }

        return new FragmentCoordReads(finalDuplicateGroups, finalSingleReads);
    }

    private static void collapseToDuplicateGroup(final DuplicateGroup duplicateGroup, final List<GroupInfo> collapsedReads)
    {
        for(GroupInfo otherGroup : collapsedReads)
        {
            if(otherGroup.Group instanceof ReadInfo)
            {
                ReadInfo otherRead = (ReadInfo) otherGroup.Group;
                duplicateGroup.addRead(otherRead.read());
            }
            else
            {
                DuplicateGroup otherDuplicateGroup = (DuplicateGroup)otherGroup.Group;
                duplicateGroup.addReads(otherDuplicateGroup.reads());
            }
        }
    }

    private void collapseGroups(final List<GroupInfo> groups)
    {
        if(groups == null || groups.size() == 1)
            return;

        Collections.sort(groups);

        while(true)
        {
            int groupIndex = findMaxGroupToCollapse(groups);

            if(groupIndex < 0)
                break;

            GroupInfo mainGroup = groups.get(groupIndex);
            List<GroupInfo> collapseGroups = findGroupsToCollapse(groups, groupIndex);

            for(GroupInfo otherGroup : collapseGroups)
            {
                mainGroup.addGroup(otherGroup);
                otherGroup.markCollapsed();
            }
        }
    }

    private int findMaxGroupToCollapse(final List<GroupInfo> groups)
    {
        int maxGroupSize = 0;
        int maxGroupIndex = 0;

        for(int i = 0; i < groups.size(); ++i)
        {
            List<GroupInfo> collapseGroups = findGroupsToCollapse(groups, i);

            if(collapseGroups == null)
                continue;

            // prioritise which would become the new largest group, not the current size
            int newSize = groups.get(i).size() + collapseGroups.stream().mapToInt(x -> x.size()).sum();

            if(newSize > maxGroupSize)
            {
                maxGroupSize = newSize;
                maxGroupIndex = i;
            }
        }

        return maxGroupSize > 0 ? maxGroupIndex : -1;
    }

    private List<GroupInfo> findGroupsToCollapse(final List<GroupInfo> groups, final int groupIndex)
    {
        GroupInfo groupInfo = groups.get(groupIndex);

        if(groupInfo.collapsed())
            return null;

        List<GroupInfo> collapseGroups = null;

        for(int i = 0; i <= 1; ++i)
        {
            boolean searchUp = (i == 0);
            int j = searchUp ? groupIndex + 1 : groupIndex - 1;

            while(j >= 0 && j < groups.size())
            {
                GroupInfo nextGroup = groups.get(j);

                if(!groupInfo.withinRange(nextGroup))
                {
                    if(nextGroup.PositionLower > groupInfo.PositionLower + mMaxDuplicateDistance)
                        break;
                }
                else if(!nextGroup.collapsed())
                {
                    if(collapseGroups == null)
                        collapseGroups = Lists.newArrayListWithCapacity(mMaxDuplicateDistance);

                    collapseGroups.add(nextGroup);
                }

                j += searchUp ? 1 : -1;
            }
        }

        return collapseGroups;
    }

    public FragmentCoordReads collapseGroupsOld(
            @Nullable final List<DuplicateGroup> duplicateGroups, @Nullable final List<ReadInfo> singleReads)
    {
        Map<String,Map<FragStartEnd,DuplicateGroup>> keyGroups = Maps.newHashMap();

        if(singleReads != null)
            singleReads.forEach(x -> addSingleRead(x, keyGroups));

        if(duplicateGroups != null)
            duplicateGroups.forEach(x -> addDuplicateGroup(x, keyGroups));

        return getCollapsedGroups(keyGroups);
    }

    private record FragStartEnd(int fragStartPos, int fragEndPos) implements Comparable<FragStartEnd>
    {
        @Override
        public int compareTo(final FragStartEnd o)
        {
            int diffFragStartPos = fragStartPos - o.fragStartPos;
            if(diffFragStartPos != 0)
            {
                return diffFragStartPos;
            }

            return fragEndPos - o.fragEndPos;
        }

        public int distance(final FragStartEnd o)
        {
            return abs(fragStartPos - o.fragStartPos) + abs(fragEndPos - o.fragEndPos);
        }
    }

    public void addSingleRead(final ReadInfo readInfo, final Map<String,Map<FragStartEnd,DuplicateGroup>> keyGroups)
    {
        DuplicateGroup duplicateGroup = new DuplicateGroup(null, readInfo.read(), readInfo.coordinates());
        addDuplicateGroup(duplicateGroup, keyGroups);
    }

    public void addDuplicateGroup(final DuplicateGroup duplicateGroup, final Map<String,Map<FragStartEnd,DuplicateGroup>> keyGroups)
    {
        FragmentCoords coords = duplicateGroup.fragmentCoordinates();
        String collapsedKey = collapseKey(coords);
        int fragStartPos = coords.ReadIsLower ? coords.PositionLower : coords.PositionUpper;
        int fragEndPos = coords.ReadIsLower ? coords.PositionUpper : coords.PositionLower;

        keyGroups.computeIfAbsent(collapsedKey, key -> Maps.newHashMap());

        keyGroups.get(collapsedKey).merge(
                new FragStartEnd(fragStartPos, fragEndPos), duplicateGroup, DUPLICATE_GROUP_MERGER);
    }

    public FragmentCoordReads getCollapsedGroups(final Map<String,Map<FragStartEnd,DuplicateGroup>> keyGroups)
    {
        if(keyGroups.isEmpty())
            return null;

        BiPredicate<FragStartEnd, FragStartEnd> canMergeFn = (x, y) -> x.distance(y) <= mMaxDuplicateDistance;
        List<DuplicateGroup> collapsedGroups = Lists.newArrayList();
        for(Map<FragStartEnd, DuplicateGroup> keyGroup : keyGroups.values())
        {
            collapsedGroups.addAll(clusterMerger(
                    keyGroup, canMergeFn, DUPLICATE_GROUP_COMPARATOR, DUPLICATE_GROUP_MERGER, null));
        }

        return CollapseUtils.getFragmentCoordReads(collapsedGroups.stream());
    }

    private static String collapseKey(final FragmentCoords fragmentCoords)
    {
        String key = fragmentCoords.OrientLower.isForward() ? "F" : "R";
        if(fragmentCoords.SuppReadInfo != null)
        {
            return key + "_S";
        }

        return key;
    }
}

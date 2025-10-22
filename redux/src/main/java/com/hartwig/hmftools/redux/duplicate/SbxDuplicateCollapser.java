package com.hartwig.hmftools.redux.duplicate;

import static java.lang.Math.abs;
import static java.lang.String.format;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.redux.common.ReadInfo;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public class SbxDuplicateCollapser
{
    private final int mMaxDuplicateDistance;

    private final Map<Integer,Integer> mCollapsingDistanceFrequency;

    public SbxDuplicateCollapser(int maxDuplicateDistance)
    {
        mMaxDuplicateDistance = maxDuplicateDistance;
        mCollapsingDistanceFrequency = Maps.newHashMap();
    }

    private static boolean withinRange(int firstPosLower, int firstPosUpper, int secondPosLower, int secondPosUpper, int maxDistance)
    {
        return abs(firstPosLower - secondPosLower) + abs(firstPosUpper - secondPosUpper) <= maxDistance;
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
            return SbxDuplicateCollapser.withinRange(
                    PositionLower, PositionUpper, other.PositionLower, other.PositionUpper, mMaxDuplicateDistance);
        }

        public int distance(final GroupInfo other)
        {
            return abs(PositionLower - other.PositionLower) + abs(PositionUpper - other.PositionUpper);
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

        public String readId()
        {
            if(Group instanceof ReadInfo)
                return ((ReadInfo)Group).id();
            else
                return ((DuplicateGroup)Group).reads().get(0).getReadName();
        }

        public String toString() { return format("coords(%d - %d) size(%d)", PositionLower, PositionUpper, mSize); }
    }

    public FragmentCoordReads collapseGroups(
            @Nullable final List<DuplicateGroup> duplicateGroups, @Nullable final List<ReadInfo> singleReads)
    {
        // add single reads and duplicate groups to collections which can be collapsed - same orientation, and primary/supplementary
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

        // run collapsing routine - this iteratively finds the group which can collapse in the most adjacent reads & other groups
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

                    // convert single read to group
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

                checkSupplementaryConsistency(duplicateGroup);

                finalDuplicateGroups.add(duplicateGroup);
            }
        }

        return new FragmentCoordReads(finalDuplicateGroups, finalSingleReads);
    }

    private static void collapseToDuplicateGroup(final DuplicateGroup duplicateGroup, final List<GroupInfo> collapsedReads)
    {
        // for SBX, collapsed reads will count towards consensus
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

            /*
            if(ReduxConfig.LogReadIds.stream().anyMatch(x -> x.equals(mainGroup.readId())))
            {
                RD_LOGGER.debug("SBX duplicate collapse: {}", mainGroup.readId());
            }
            */

            List<GroupInfo> collapseGroups = findGroupsToCollapse(groups, groupIndex);

            for(GroupInfo otherGroup : collapseGroups)
            {
                mainGroup.addGroup(otherGroup);
                otherGroup.markCollapsed();

                /*
                int collapseDistance = mainGroup.distance(otherGroup);
                addStats(collapseDistance);
                */
            }
        }
    }

    private void addStats(int collapseDistance)
    {
        mCollapsingDistanceFrequency.put(collapseDistance, mCollapsingDistanceFrequency.getOrDefault(collapseDistance, 0) + 1);
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

    private static String collapseKey(final FragmentCoords fragmentCoords)
    {
        String key = fragmentCoords.OrientLower.isForward() ? "F" : "R";
        if(fragmentCoords.SuppReadInfo != null)
        {
            return key + "_S";
        }

        return key;
    }

    private void checkSupplementaryConsistency(final DuplicateGroup duplicateGroup)
    {
        // supplementaries are not guaranteed to align, especially if the primaries were grouped with jitter
        // in which case discard any read which is outside the valid jitter range to help with consensus

        List<SAMRecord> suppReads = duplicateGroup.reads().stream().filter(x -> x.getSupplementaryAlignmentFlag()).collect(Collectors.toList());

        if(suppReads.isEmpty())
            return;

        SAMRecord templateRead;

        if(suppReads.size() == duplicateGroup.reads().size())
        {
            // could instead choose the most common by cigar or frag coords
            templateRead = duplicateGroup.reads().get(0);
        }
        else
        {
            templateRead = duplicateGroup.reads().stream().filter(x -> !x.getSupplementaryAlignmentFlag()).findFirst().orElse(null);

            if(templateRead == null)
                templateRead = duplicateGroup.reads().get(0);
        }

        // discard from consensus any read not within range of the template read's coords
        FragmentCoords templateCoords = FragmentCoords.fromRead(templateRead, true, false);

        List<SAMRecord> nonConsensusReads = null;

        for(SAMRecord read : duplicateGroup.reads())
        {
            if(read == templateRead)
                continue;

            FragmentCoords readCoords = FragmentCoords.fromRead(read, true, false);

            if(!SbxDuplicateCollapser.withinRange(
                    templateCoords.PositionLower, templateCoords.PositionUpper, readCoords.PositionLower, readCoords.PositionUpper,
                    mMaxDuplicateDistance))
            {
                if(nonConsensusReads == null)
                    nonConsensusReads = Lists.newArrayList();

                nonConsensusReads.add(read);
            }
        }

        if(nonConsensusReads != null)
        {
            nonConsensusReads.forEach(x -> duplicateGroup.reads().remove(x));
            duplicateGroup.addNonConsensusReads(nonConsensusReads);
        }
    }
}

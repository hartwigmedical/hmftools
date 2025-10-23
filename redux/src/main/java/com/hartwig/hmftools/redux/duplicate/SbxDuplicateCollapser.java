package com.hartwig.hmftools.redux.duplicate;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.util.Collections;
import java.util.Comparator;
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

    private static final int MAX_CHAINED_MULTIPLE = 3;

    private final Map<Integer,Integer> mCollapsingDistanceFrequency;

    public SbxDuplicateCollapser(int maxDuplicateDistance)
    {
        mMaxDuplicateDistance = maxDuplicateDistance;
        mCollapsingDistanceFrequency = Maps.newHashMap();
    }

    /* collapsing routine and conditions
        - form a group of single and duplicate groups read for coordinate comparison, ordered by lower then upper unclipped positions
        - primary vs supplementary and forward vs reverse orientation reads are collapsed separately
        - find all non-overlapping groups to collapse in proximate smaller groups
        - process these, taking the largest first
        - take the wider coordinates from collapsed groups, ie allow the lower and upper positions to expand to allow further collapsing
        - continue until all collapsing is exhausted

    */

    private static boolean withinRange(int firstPosLower, int firstPosUpper, int secondPosLower, int secondPosUpper, int maxDistance)
    {
        return abs(firstPosLower - secondPosLower) + abs(firstPosUpper - secondPosUpper) <= maxDistance;
    }

    private class GroupInfo implements Comparable<GroupInfo>
    {
        public final Object Group;

        public final int PositionLower;
        public final int PositionUpper;

        public int mPositionLowerMin;
        public int mPositionLowerMax;
        private int mPositionUpperMin;
        private int mPositionUpperMax;

        private int mIndex; // in its group list, and is immutable
        private int mSize;
        private boolean mCollapsed;
        private List<GroupInfo> mCollapsedGroups;

        public GroupInfo(final Object group, final int positionLower, final int positionUpper, final int initialSize)
        {
            mIndex = -1;
            Group = group;
            PositionLower = positionLower;
            PositionUpper = positionUpper;

            mPositionLowerMin = positionLower;
            mPositionLowerMax = positionLower;

            mPositionUpperMin = positionUpper;
            mPositionUpperMax = positionUpper;

            mSize = initialSize;
            mCollapsedGroups = null;
            mCollapsed = false;
        }

        public void setIndex(int index) { mIndex = index; }
        public int index() { return mIndex; }

        public int posLowerMin() { return mPositionLowerMin; }
        public int posLowerMax() { return mPositionLowerMax; }
        public int posUpperMin() { return mPositionUpperMin; }
        public int posUpperMax() { return mPositionUpperMax; }

        public int size() { return mSize; }
        public void markCollapsed() { mCollapsed = true; }
        public boolean collapsed() { return mCollapsed; }

        public boolean hasCollapsedGroups() { return mCollapsedGroups != null; }
        public List<GroupInfo> collapsedGroups() { return mCollapsedGroups; }

        public void addGroup(final GroupInfo group)
        {
            if(mCollapsedGroups == null)
                mCollapsedGroups = Lists.newArrayListWithCapacity(2);

            mCollapsedGroups.add(group);
            mSize += group.size();

            mPositionLowerMin = min(mPositionLowerMin, group.posLowerMin());
            mPositionLowerMax = max(mPositionLowerMax, group.posLowerMax());

            mPositionUpperMin = min(mPositionUpperMin, group.posUpperMin());
            mPositionUpperMax = max(mPositionUpperMax, group.posUpperMax());

            // cap how far the chained collapsing can extend
            int maxDistanceFromOriginal = MAX_CHAINED_MULTIPLE * mMaxDuplicateDistance;
            mPositionLowerMin = max(mPositionLowerMin, PositionLower - maxDistanceFromOriginal);
            mPositionLowerMax = min(mPositionLowerMax, PositionLower + maxDistanceFromOriginal);
            mPositionUpperMin = max(mPositionUpperMin, PositionUpper - maxDistanceFromOriginal);
            mPositionUpperMax = min(mPositionUpperMax, PositionUpper + maxDistanceFromOriginal);
        }

        public List<SAMRecord> allReads()
        {
            // gathers reads including recursively from collapsed groups
            List<SAMRecord> reads = Lists.newArrayListWithCapacity(mSize);

            if(Group instanceof ReadInfo)
            {
                ReadInfo readInfo = (ReadInfo)Group;
                reads.add(readInfo.read());
            }
            else
            {
                DuplicateGroup duplicateGroup = (DuplicateGroup)Group;
                reads.addAll(duplicateGroup.reads());
            }

            if(mCollapsedGroups != null)
            {
                for(GroupInfo collapsedGroup : mCollapsedGroups)
                {
                    reads.addAll(collapsedGroup.allReads());
                }
            }

            return reads;
        }

        private static int minDistance(int firstMin, int firstMax, int secondMin, int secondMax)
        {
            if(firstMax < secondMin)
                return secondMin - firstMax;
            else if(secondMax < firstMin)
                return firstMin - secondMax;
            else
                return 0;
        }

        public boolean withinRange(final GroupInfo other)
        {
            return distance(other) <= mMaxDuplicateDistance;
        }

        public int distance(final GroupInfo other)
        {
            int lowerDistance = minDistance(mPositionLowerMin, mPositionLowerMax, other.posLowerMin(), other.posLowerMax());
            int upperDistance = minDistance(mPositionUpperMin, mPositionUpperMax, other.posUpperMin(), other.posUpperMax());
            return lowerDistance + upperDistance;
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

        public String toString()
        {
            StringBuilder sb = new StringBuilder(format("%d: coords(%d-%d - %d-%d)",
                    mIndex, mPositionLowerMin, mPositionLowerMax, mPositionUpperMin, mPositionUpperMax));

            if(mSize > 0)
            {
                sb.append(format(" size(%d)", mSize));

                if(mCollapsedGroups != null)
                    sb.append(format(" collapsedGroups(%d)", mCollapsedGroups.size()));
            }
            else
            {
                sb.append(" collapsed");
            }

            return sb.toString();
        }
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
            if(groups.size() > 1)
            {
                Collections.sort(groups); // by start then end position

                for(int i = 0; i < groups.size(); ++i)
                {
                    groups.get(i).setIndex(i);
                }

                collapseGroups(groups);
            }

            for(GroupInfo groupInfo : groups)
            {
                if(groupInfo.collapsed())
                    continue;

                /*
                if(groupInfo.posLowerMax() - groupInfo.posLowerMin() > mMaxDuplicateDistance * 20
                || groupInfo.posUpperMax() - groupInfo.posUpperMin() > mMaxDuplicateDistance * 20)
                {
                    RD_LOGGER.debug("groupInfo({}) has wider position range", groupInfo);
                }
                */

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

    private static void collapseToDuplicateGroup(final DuplicateGroup duplicateGroup, final List<GroupInfo> collapsedGroups)
    {
        // for SBX, collapsed reads will count towards consensus
        for(GroupInfo otherGroup : collapsedGroups)
        {
            duplicateGroup.addReads(otherGroup.allReads());
        }
    }

    private void collapseGroups(final List<GroupInfo> groups)
    {
        while(true)
        {
            List<GroupInfo> groupsToCollapse = findGroupsToCollapse(groups);

            if(groupsToCollapse == null || groupsToCollapse.isEmpty())
                break;

            // process in order of largest groups first
            Collections.sort(groupsToCollapse, Comparator.comparingInt(x -> -x.size()));

            for(GroupInfo groupInfo : groupsToCollapse)
            {
                if(groupInfo.collapsed()) // may have just been collapsed
                    continue;

                /*
                if(ReduxConfig.LogReadIds.stream().anyMatch(x -> x.equals(mainGroup.readId())))
                {
                    RD_LOGGER.debug("SBX duplicate collapse: {}", mainGroup.readId());
                }
                */

                List<GroupInfo> collapseGroups = findProximateGroups(groups, groupInfo);

                if(collapseGroups == null || collapseGroups.isEmpty())
                    continue;

                for(GroupInfo otherGroup : collapseGroups)
                {
                    if(otherGroup.collapsed())
                        continue;

                    groupInfo.addGroup(otherGroup);
                    otherGroup.markCollapsed();

                    /*
                    int collapseDistance = mainGroup.distance(otherGroup);
                    addStats(collapseDistance);
                    */
                }
            }
        }
    }

    private List<GroupInfo> findGroupsToCollapse(final List<GroupInfo> groups)
    {
        GroupInfo lastCollapsingGroup = null;

        List<GroupInfo> groupToCollapse = null;

        for(int i = 0; i < groups.size(); ++i)
        {
            GroupInfo groupInfo = groups.get(i);

            if(groupInfo.collapsed())
                continue;

            if(lastCollapsingGroup != null)
            {
                // avoid checking a group which is close to the previous groups which collapsed in others unless it's larger
                if(groupInfo.size() <= lastCollapsingGroup.size()
                && withinRange(
                        lastCollapsingGroup.PositionLower, lastCollapsingGroup.PositionUpper,
                        groupInfo.PositionLower, groupInfo.PositionUpper, mMaxDuplicateDistance * 2))
                {
                    continue;
                }
            }

            List<GroupInfo> collapseGroups = findProximateGroups(groups, groupInfo);

            if(collapseGroups == null)
                continue;

            // take this group's proximate collapsed groups if it isn't within range of another larger group
            if(groupToCollapse == null)
                groupToCollapse = Lists.newArrayList();

            groupToCollapse.add(groupInfo);
            lastCollapsingGroup = groupInfo;
        }

        return groupToCollapse;
    }

    private List<GroupInfo> findProximateGroups(final List<GroupInfo> groups, final GroupInfo groupInfo)
    {
        List<GroupInfo> collapseGroups = null;

        for(int i = 0; i <= 1; ++i)
        {
            boolean searchUp = (i == 0);
            int j = searchUp ? groupInfo.index() + 1 : groupInfo.index() - 1;

            while(j >= 0 && j < groups.size())
            {
                GroupInfo nextGroup = groups.get(j);

                if(!groupInfo.withinRange(nextGroup))
                {
                    if(nextGroup.PositionLower > groupInfo.PositionLower + mMaxDuplicateDistance)
                        break;
                }
                else if(!nextGroup.collapsed() && nextGroup.size() <= groupInfo.size())
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

    private void addStats(int collapseDistance)
    {
        mCollapsingDistanceFrequency.put(collapseDistance, mCollapsingDistanceFrequency.getOrDefault(collapseDistance, 0) + 1);
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

package com.hartwig.hmftools.redux.duplicate;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.collect.MergeUtils.clusterMerger;
import static com.hartwig.hmftools.redux.duplicate.CollapseUtils.DUPLICATE_GROUP_MERGER;
import static com.hartwig.hmftools.redux.duplicate.DuplicateGroup.DUPLICATE_GROUP_COMPARATOR;

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

    public FragmentCoordReads collapseGroups(
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
        String key = fragmentCoords.ReadIsLower ? "" : "R";
        if(fragmentCoords.SuppReadInfo != null)
        {
            return key + "S";
        }

        return key;
    }
}

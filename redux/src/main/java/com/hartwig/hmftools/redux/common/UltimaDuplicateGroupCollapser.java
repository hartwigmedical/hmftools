package com.hartwig.hmftools.redux.common;

import java.util.List;
import java.util.Map;
import java.util.SortedMap;

import com.google.common.collect.Maps;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public class UltimaDuplicateGroupCollapser extends DuplicateGroupCollapser
{
    public static int MAX_THREE_PRIME_COLLAPSE_DISTANCE = 10;

    private static class GroupCollapser
    {
        private final Map<String, SortedMap<Integer, List<ReadWithFragCoords>>> mGroups;

        public GroupCollapser()
        {
            mGroups = Maps.newHashMap();
        }

        public void addRead(final SAMRecord read, final FragmentCoords fragmentCoords)
        {
            List<ReadWithFragCoords> readGroup = getOrCreateReadGroup(fragmentCoords);
            readGroup.add(new ReadWithFragCoords(read, fragmentCoords));
        }

        public FragmentCoordReads getCollapsedGroups()
        {
            if(mGroups.isEmpty())
                return null;

            List<ReadInfo> singleReads = Lists.newArrayList();
            List<DuplicateGroup> duplicateGroups = Lists.newArrayList();

            for(SortedMap<Integer, List<ReadWithFragCoords>> fivePrimeGroup : mGroups.values())
            {
                List<ReadWithFragCoords> collapedGroup = null;
                int lastFragEndPos = 0;
                for(Map.Entry<Integer, List<ReadWithFragCoords>> fivePrimeGroupEntry : fivePrimeGroup.entrySet())
                {
                    int fragEndPos = fivePrimeGroupEntry.getKey();
                    List<ReadWithFragCoords> reads = fivePrimeGroupEntry.getValue();

                    if(collapedGroup == null)
                    {
                        lastFragEndPos = fragEndPos;
                        collapedGroup = Lists.newArrayList();
                        collapedGroup.addAll(reads);
                        continue;
                    }

                    if(fragEndPos - lastFragEndPos <= MAX_THREE_PRIME_COLLAPSE_DISTANCE)
                    {
                        lastFragEndPos = fragEndPos;
                        collapedGroup.addAll(reads);
                        continue;
                    }

                    if(collapedGroup.size() == 1)
                    {
                        ReadWithFragCoords read = collapedGroup.get(0);
                        singleReads.add(new ReadInfo(read));
                    }
                    else
                    {
                        duplicateGroups.add(new MultiCoordDuplicateGroup(collapedGroup));
                    }

                    lastFragEndPos = fragEndPos;
                    collapedGroup = Lists.newArrayList();
                    collapedGroup.addAll(reads);
                }

                if(collapedGroup == null && collapedGroup.isEmpty())
                    continue;

                if(collapedGroup.size() == 1)
                {
                    ReadWithFragCoords read = collapedGroup.get(0);
                    singleReads.add(new ReadInfo(read));
                }
                else
                {
                    duplicateGroups.add(new MultiCoordDuplicateGroup(collapedGroup));
                }
            }

            return new FragmentCoordReads(duplicateGroups, singleReads);
        }

        private SortedMap<Integer, List<ReadWithFragCoords>> getOrCreateFivePrimeGroup(final String fivePrimeKey)
        {
            SortedMap<Integer, List<ReadWithFragCoords>> fivePrimeGroup = mGroups.get(fivePrimeKey);
            if(fivePrimeGroup != null)
                return fivePrimeGroup;

            fivePrimeGroup = Maps.newTreeMap();
            mGroups.put(fivePrimeKey, fivePrimeGroup);
            return fivePrimeGroup;
        }

        private List<ReadWithFragCoords> getOrCreateReadGroup(final FragmentCoords fragmentCoords)
        {
            String fivePrimeKey = collapseToFivePrimeKey(fragmentCoords);
            SortedMap<Integer, List<ReadWithFragCoords>> fivePrimeGroup = getOrCreateFivePrimeGroup(fivePrimeKey);
            int fragEndPos = fragmentCoords.ReadIsLower ? fragmentCoords.PositionUpper : fragmentCoords.PositionLower;
            List<ReadWithFragCoords> readGroup = fivePrimeGroup.get(fragEndPos);
            if(readGroup != null)
                return readGroup;

            readGroup = Lists.newArrayList();
            fivePrimeGroup.put(fragEndPos, readGroup);
            return readGroup;
        }
    }

    public UltimaDuplicateGroupCollapser()
    {
    }

    @Override
    public FragmentCoordReads collapse(@Nullable final List<DuplicateGroup> duplicateGroups, @Nullable final List<ReadInfo> singleReads)
    {
        GroupCollapser groupCollapser = new GroupCollapser();

        if(singleReads != null)
        {
            for(ReadInfo readInfo : singleReads)
                groupCollapser.addRead(readInfo.read(), readInfo.coordinates());
        }

        if(duplicateGroups != null)
        {
            for(DuplicateGroup duplicateGroup : duplicateGroups)
                duplicateGroup.processReads(groupCollapser::addRead);
        }

        return groupCollapser.getCollapsedGroups();
    }

    private static String collapseToFivePrimeKey(final FragmentCoords fragmentCoords)
    {
        String key;
        if(fragmentCoords.ReadIsLower)
        {
            key = fragmentCoords.ChromsomeLower + ":" + fragmentCoords.PositionLower;
        }
        else
        {
            key = fragmentCoords.ChromsomeUpper + ":" + fragmentCoords.PositionUpper + "_R";
        }

        if(fragmentCoords.SuppReadInfo != null)
            return key + "_S";

        return key;
    }
}

package com.hartwig.hmftools.redux.common;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.sequencing.SequencingType.BIOMODAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ULTIMA;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.BinaryOperator;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.collect.Heap;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

@FunctionalInterface
public interface DuplicateGroupCollapser
{
    FragmentCoordReads collapse(@Nullable final List<DuplicateGroup> duplicateGroups, @Nullable final List<ReadInfo> singleReads);

    static DuplicateGroupCollapser from(final DuplicateGroupCollapseConfig config)
    {
        if(config.Sequencing == ULTIMA)
            return DuplicateGroupCollapser::ultimaCollapse;

        if(config.Sequencing == BIOMODAL)
            return DuplicateGroupCollapser::biomodalCollapse;

        if(config.Sequencing == SBX && config.SbxMaxDuplicateDistance > 0)
            return sbxCollapserFactory(config.SbxMaxDuplicateDistance);

        return null;
    }

    static boolean isEnabled(final DuplicateGroupCollapseConfig config)
    {
        if(config.Sequencing == ULTIMA)
            return true;

        if(config.Sequencing == BIOMODAL)
            return true;

        if(config.Sequencing == SBX && config.SbxMaxDuplicateDistance > 0)
            return true;

        return false;
    }

    BinaryOperator<DuplicateGroup> DUPLICATE_GROUP_MERGER = (acc, group) ->
    {
        acc.addReads(group.reads());
        return acc;
    };

    private static FragmentCoordReads getFragmentCoordReads(final Collection<DuplicateGroup> collapsedGroups)
    {
        List<DuplicateGroup> duplicateGroups = Lists.newArrayList();
        List<ReadInfo> singleReads = Lists.newArrayList();
        for(DuplicateGroup collapsedGroup : collapsedGroups)
        {
            if(collapsedGroup.readCount() == 1)
            {
                ReadInfo singleRead = new ReadInfo(collapsedGroup.reads().get(0), collapsedGroup.fragmentCoordinates());
                singleReads.add(singleRead);
                continue;
            }

            duplicateGroups.add(collapsedGroup);
        }

        return new FragmentCoordReads(duplicateGroups, singleReads);
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

    int SINGLE_END_JITTER_COLLAPSE_DISTANCE = 10;

    class UltimaCollapser
    {
        private static final Comparator<Map.Entry<Integer, DuplicateGroup>> FIVE_PRIME_GROUP_ENTRY_COMPARATOR =
                Comparator.comparingInt((Map.Entry<Integer, DuplicateGroup> x) -> x.getValue().readCount()).reversed();

        private final Map<String, TreeMap<Integer, DuplicateGroup>> mFivePrimeGroups;

        public UltimaCollapser()
        {
            mFivePrimeGroups = Maps.newHashMap();
        }

        public void addSingleRead(final ReadInfo readInfo)
        {
            addDuplicateGroup(new DuplicateGroup(null, readInfo.read(), readInfo.coordinates()));
        }

        public void addDuplicateGroup(final DuplicateGroup duplicateGroup)
        {
            FragmentCoords coords = duplicateGroup.fragmentCoordinates();
            int fragEndPos = coords.ReadIsLower ? coords.PositionUpper : coords.PositionLower;
            String fivePrimeKey = collapseToFivePrimeKey(coords);
            TreeMap<Integer, DuplicateGroup> fivePrimeGroup = getOrCreateFivePrimeGroup(fivePrimeKey);
            fivePrimeGroup.merge(fragEndPos, duplicateGroup, DUPLICATE_GROUP_MERGER);
        }

        public FragmentCoordReads getCollapsedGroups()
        {
            if(mFivePrimeGroups.isEmpty())
                return null;

            List<DuplicateGroup> collapsedGroups = Lists.newArrayList();
            for(TreeMap<Integer, DuplicateGroup> fivePrimeGroup : mFivePrimeGroups.values())
            {
                Heap<Map.Entry<Integer, DuplicateGroup>> dupGroupHeap = new Heap<>(FIVE_PRIME_GROUP_ENTRY_COMPARATOR);
                dupGroupHeap.addAll(fivePrimeGroup.entrySet());
                while(!dupGroupHeap.isEmpty())
                {
                    Map.Entry<Integer, DuplicateGroup> entry = dupGroupHeap.pop();
                    if(!fivePrimeGroup.containsKey(entry.getKey()))
                        continue;

                    int baseKey = entry.getKey();
                    int minKey = baseKey;
                    int maxKey = baseKey;
                    DuplicateGroup collapsedGroup = entry.getValue();
                    fivePrimeGroup.remove(baseKey);
                    NavigableMap<Integer, DuplicateGroup> beforeBase = fivePrimeGroup.headMap(baseKey, false).descendingMap();
                    NavigableMap<Integer, DuplicateGroup> afterBase = fivePrimeGroup.tailMap(baseKey, false);
                    while(!beforeBase.isEmpty() || !afterBase.isEmpty())
                    {
                        Map.Entry<Integer, DuplicateGroup> beforeEntry = beforeBase.firstEntry();
                        Map.Entry<Integer, DuplicateGroup> afterEntry = afterBase.firstEntry();
                        if(beforeEntry == null)
                        {
                            int distance = afterEntry.getKey() - minKey;
                            if(distance > SINGLE_END_JITTER_COLLAPSE_DISTANCE)
                                break;

                            maxKey = afterEntry.getKey();
                            collapsedGroup.addReads(afterEntry.getValue().reads());
                            afterBase.pollFirstEntry();
                            continue;
                        }

                        if(afterEntry == null)
                        {
                            int distance = maxKey - beforeEntry.getKey();
                            if(distance > SINGLE_END_JITTER_COLLAPSE_DISTANCE)
                                break;

                            minKey = beforeEntry.getKey();
                            collapsedGroup.addReads(beforeEntry.getValue().reads());
                            beforeBase.pollFirstEntry();
                            continue;
                        }

                        int afterDistance = afterEntry.getKey() - minKey;
                        int beforeDistance = maxKey - beforeEntry.getKey();
                        if(afterDistance > SINGLE_END_JITTER_COLLAPSE_DISTANCE && beforeDistance > SINGLE_END_JITTER_COLLAPSE_DISTANCE)
                            break;

                        boolean useBeforeEntry = beforeDistance <= SINGLE_END_JITTER_COLLAPSE_DISTANCE && (
                                afterDistance > SINGLE_END_JITTER_COLLAPSE_DISTANCE
                                        || beforeEntry.getValue().readCount() >= afterEntry.getValue().readCount());
                        if(useBeforeEntry)
                        {
                            minKey = beforeEntry.getKey();
                            collapsedGroup.addReads(beforeEntry.getValue().reads());
                            beforeBase.pollFirstEntry();
                            continue;
                        }

                        maxKey = afterEntry.getKey();
                        collapsedGroup.addReads(afterEntry.getValue().reads());
                        afterBase.pollFirstEntry();
                    }

                    collapsedGroups.add(collapsedGroup);
                }
            }

            return getFragmentCoordReads(collapsedGroups);
        }

        private TreeMap<Integer, DuplicateGroup> getOrCreateFivePrimeGroup(final String fivePrimeKey)
        {
            mFivePrimeGroups.computeIfAbsent(fivePrimeKey, key -> Maps.newTreeMap());
            return mFivePrimeGroups.get(fivePrimeKey);
        }
    }

    static FragmentCoordReads ultimaCollapse(
            @Nullable final List<DuplicateGroup> duplicateGroups, @Nullable final List<ReadInfo> singleReads)
    {
        UltimaCollapser collapser = new UltimaCollapser();

        if(singleReads != null)
            singleReads.forEach(collapser::addSingleRead);

        if(duplicateGroups != null)
            duplicateGroups.forEach(collapser::addDuplicateGroup);

        return collapser.getCollapsedGroups();
    }

    static FragmentCoordReads biomodalCollapse(
            @Nullable final List<DuplicateGroup> duplicateGroups, @Nullable final List<ReadInfo> singleReads)
    {
        Map<String, DuplicateGroup> fivePrimeGroups = Maps.newHashMap();

        if(singleReads != null)
        {
            for(ReadInfo readInfo : singleReads)
            {
                SAMRecord read = readInfo.read();
                FragmentCoords coords = readInfo.coordinates();
                String fivePrimeKey = collapseToFivePrimeKey(coords);
                DuplicateGroup duplicateGroup = new DuplicateGroup(null, read, coords);
                fivePrimeGroups.merge(fivePrimeKey, duplicateGroup, DUPLICATE_GROUP_MERGER);
            }
        }

        if(duplicateGroups != null)
        {
            for(DuplicateGroup duplicateGroup : duplicateGroups)
            {
                String fivePrimeKey = collapseToFivePrimeKey(duplicateGroup.fragmentCoordinates());
                fivePrimeGroups.merge(fivePrimeKey, duplicateGroup, DUPLICATE_GROUP_MERGER);
            }
        }

        if(fivePrimeGroups.isEmpty())
            return null;

        return getFragmentCoordReads(fivePrimeGroups.values());
    }

    class SbxCollapser
    {
        private record FragStartEnd(int fragStartPos, int fragEndPos) implements Comparable<FragStartEnd>
        {
            @Override
            public int compareTo(final FragStartEnd o)
            {
                int diffFragStartPos = fragStartPos - o.fragStartPos;
                if(diffFragStartPos != 0)
                    return diffFragStartPos;

                return fragEndPos - o.fragEndPos;
            }

            public int distance(final FragStartEnd o)
            {
                return abs(fragStartPos - o.fragStartPos) + abs(fragEndPos - o.fragEndPos);
            }
        }

        private static final Comparator<Map.Entry<FragStartEnd, DuplicateGroup>> KEY_GROUP_ENTRY_COMPARATOR =
                Comparator.comparingInt((Map.Entry<FragStartEnd, DuplicateGroup> x) -> x.getValue().readCount()).reversed();

        private final int mMaxDuplicateDistance;
        private final Map<String, Map<FragStartEnd, DuplicateGroup>> mKeyGroups;

        public SbxCollapser(int maxDuplicateDistance)
        {
            mMaxDuplicateDistance = maxDuplicateDistance;
            mKeyGroups = Maps.newHashMap();
        }

        public void addSingleRead(final ReadInfo readInfo)
        {
            DuplicateGroup duplicateGroup = new DuplicateGroup(null, readInfo.read(), readInfo.coordinates());
            addDuplicateGroup(duplicateGroup);
        }

        public void addDuplicateGroup(final DuplicateGroup duplicateGroup)
        {
            FragmentCoords coords = duplicateGroup.fragmentCoordinates();
            String collapsedKey = collapseKey(coords);
            int fragStartPos = coords.ReadIsLower ? coords.PositionLower : coords.PositionUpper;
            int fragEndPos = coords.ReadIsLower ? coords.PositionUpper : coords.PositionLower;
            mKeyGroups.computeIfAbsent(collapsedKey, key -> Maps.newHashMap());
            mKeyGroups.get(collapsedKey).merge(new FragStartEnd(fragStartPos, fragEndPos), duplicateGroup, DUPLICATE_GROUP_MERGER);
        }

        public FragmentCoordReads getCollapsedGroups()
        {
            if(mKeyGroups.isEmpty())
                return null;

            List<DuplicateGroup> collapsedGroups = Lists.newArrayList();
            for(Map<FragStartEnd, DuplicateGroup> keyGroup : mKeyGroups.values())
            {
                if(keyGroup.size() == 1)
                {
                    collapsedGroups.addAll(keyGroup.values());
                    continue;
                }

                List<FragStartEnd> allKeys = Lists.newArrayList(keyGroup.keySet());
                Map<FragStartEnd, Set<FragStartEnd>> keyAdjacency = Maps.newHashMap();
                for(FragStartEnd key : allKeys)
                    keyAdjacency.put(key, Sets.newHashSet());

                for(int i = 0; i < allKeys.size() - 1; i++)
                {
                    FragStartEnd key1 = allKeys.get(i);
                    for(int j = i + 1; j < allKeys.size(); j++)
                    {
                        FragStartEnd key2 = allKeys.get(j);
                        if(key1.distance(key2) <= mMaxDuplicateDistance)
                        {
                            keyAdjacency.get(key1).add(key2);
                            keyAdjacency.get(key2).add(key1);
                        }
                    }
                }

                Heap<Map.Entry<FragStartEnd, DuplicateGroup>> entryHeap = new Heap<>(KEY_GROUP_ENTRY_COMPARATOR);
                entryHeap.addAll(keyGroup.entrySet());
                while(!entryHeap.isEmpty())
                {
                    Map.Entry<FragStartEnd, DuplicateGroup> entry = entryHeap.pop();
                    if(!keyGroup.containsKey(entry.getKey()))
                        continue;

                    keyGroup.remove(entry.getKey());
                    DuplicateGroup collapsedGroup = entry.getValue();
                    Set<FragStartEnd> neighbours = keyAdjacency.get(entry.getKey()).stream()
                            .filter(keyGroup::containsKey)
                            .collect(Collectors.toCollection(Sets::newHashSet));

                    if(neighbours.isEmpty())
                    {
                        collapsedGroups.add(collapsedGroup);
                        continue;
                    }

                    Heap<Map.Entry<FragStartEnd, DuplicateGroup>> neighbourHeap = new Heap<>(KEY_GROUP_ENTRY_COMPARATOR);
                    for(FragStartEnd neighbour : neighbours)
                    {
                        Map.Entry<FragStartEnd, DuplicateGroup> neighbourEntry = Pair.of(neighbour, keyGroup.get(neighbour));
                        neighbourHeap.add(neighbourEntry);
                    }

                    while(!neighbourHeap.isEmpty())
                    {
                        Map.Entry<FragStartEnd, DuplicateGroup> neighbourEntry = neighbourHeap.pop();
                        if(!neighbours.contains(neighbourEntry.getKey()))
                            continue;

                        keyGroup.remove(neighbourEntry.getKey());
                        neighbours.retainAll(keyAdjacency.get(neighbourEntry.getKey()));
                        collapsedGroup.addReads(neighbourEntry.getValue().reads());
                    }

                    collapsedGroups.add(collapsedGroup);
                }
            }

            return getFragmentCoordReads(collapsedGroups);
        }

        private static String collapseKey(final FragmentCoords fragmentCoords)
        {
            String key = fragmentCoords.ReadIsLower ? "" : "R";
            if(fragmentCoords.SuppReadInfo != null)
                return key + "S";

            return key;
        }
    }

    static DuplicateGroupCollapser sbxCollapserFactory(final int maxDuplicateDistance)
    {
        return (duplicateGroups, singleReads) ->
        {
            SbxCollapser collapser = new SbxCollapser(maxDuplicateDistance);

            if(singleReads != null)
                singleReads.forEach(collapser::addSingleRead);

            if(duplicateGroups != null)
                duplicateGroups.forEach(collapser::addDuplicateGroup);

            return collapser.getCollapsedGroups();
        };
    }
}

package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.sequencing.SequencingType.BIOMODAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ILLUMINA;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ULTIMA;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.function.BinaryOperator;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.collect.OneDGridMap;
import com.hartwig.hmftools.common.collect.TwoDGridMap;
import com.hartwig.hmftools.common.collect.UnionFind;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

@FunctionalInterface
public interface DuplicateGroupCollapser
{
    FragmentCoordReads collapse(@Nullable final List<DuplicateGroup> duplicateGroups, @Nullable final List<ReadInfo> singleReads);

    static DuplicateGroupCollapser from(final DuplicateGroupCollapseConfig config)
    {
        if(config.Sequencing == ILLUMINA && config.UMIEnabled)
            return DuplicateGroupCollapser::illuminaCollapse;

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
        if(config.Sequencing == ILLUMINA && config.UMIEnabled)
            return true;

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

    @VisibleForTesting
    @Nullable
    static String collapseToKeyWithoutCoordinates(final FragmentCoords fragmentCoords, boolean keyByFragmentOrientation)
    {
        if(fragmentCoords.Unpaired)
            return null;

        String lowerOrientation = fragmentCoords.OrientLower == FORWARD ? "F" : "R";
        String suppSuffix = fragmentCoords.SuppReadInfo == null ? "" : ":S";
        if(fragmentCoords.PositionUpper == NO_POSITION)
        {
            String unmappedSuffix = fragmentCoords.UnmappedSourced ? ":U" : "";
            return format("%s:%s%s%s", fragmentCoords.ChromsomeLower, lowerOrientation, unmappedSuffix, suppSuffix);
        }

        String upperOrientation = fragmentCoords.OrientUpper == FORWARD ? "F" : "R";
        String isLowerString = fragmentCoords.ReadIsLower ? "L" : "U";
        String fragmentOrientationSuffix = !keyByFragmentOrientation || fragmentCoords.forwardFragment() ? "" : ":N";
        return format("%s:%s:%s:%s:%s%s%s", fragmentCoords.ChromsomeLower, lowerOrientation, fragmentCoords.ChromsomeUpper, upperOrientation, isLowerString, suppSuffix, fragmentOrientationSuffix);
    }

    int SINGLE_END_JITTER_COLLAPSE_DISTANCE = 10;

    class IlluminaCollapser
    {
        private final List<DuplicateGroup> mFullyMappedGroups;
        private final Map<String, SortedMap<Integer, OneDGridMap<Integer>>> mFixedLowerGroupIndices;
        private final Map<String, SortedMap<Integer, OneDGridMap<Integer>>> mFixedUpperGroupIndices;
        private final Map<String, OneDGridMap<DuplicateGroup>> mGroupsWithUnmappedReadOrMate;

        private final Map<String, SortedMap<Integer, OneDGridMap<Integer>>> mFixedLowerCoordIndices;
        private final Map<String, SortedMap<Integer, OneDGridMap<Integer>>> mFixedUpperCoordIndices;
        private final Map<String, OneDGridMap<FragmentCoords>> mCoordsWithUnmappedReadOrMate;

        private final UnionFind<FragmentCoords> mDuplexCoordsMerger;

        public IlluminaCollapser()
        {
            mFullyMappedGroups = Lists.newArrayList();
            mFixedLowerGroupIndices = Maps.newHashMap();
            mFixedUpperGroupIndices = Maps.newHashMap();
            mGroupsWithUnmappedReadOrMate = Maps.newHashMap();

            mFixedLowerCoordIndices = Maps.newHashMap();
            mFixedUpperCoordIndices = Maps.newHashMap();
            mCoordsWithUnmappedReadOrMate = Maps.newHashMap();

            mDuplexCoordsMerger = new UnionFind<>();
        }

        public void addSingleRead(final ReadInfo readInfo)
        {
            addDuplicateGroup(new DuplicateGroup(null, readInfo.read(), readInfo.coordinates()));
        }

        public void addDuplicateGroup(final DuplicateGroup duplicateGroup)
        {
            FragmentCoords coords = duplicateGroup.fragmentCoordinates();
            mDuplexCoordsMerger.add(coords);
            String collapsedKey = collapseToKeyWithoutCoordinates(coords, true);
            String collapsedDuplexKey = collapseToKeyWithoutCoordinates(coords, false);
            if(coords.PositionUpper == NO_POSITION)
            {
                mGroupsWithUnmappedReadOrMate.computeIfAbsent(collapsedKey, key -> new OneDGridMap<>());
                mGroupsWithUnmappedReadOrMate.get(collapsedKey).put(coords.PositionLower, duplicateGroup);

                mCoordsWithUnmappedReadOrMate.computeIfAbsent(collapsedDuplexKey, key -> new OneDGridMap<>());
                mCoordsWithUnmappedReadOrMate.get(collapsedDuplexKey).put(coords.PositionLower, coords);
                return;
            }

            mFullyMappedGroups.add(duplicateGroup);
            int groupIndex = mFullyMappedGroups.size() - 1;

            mFixedLowerGroupIndices.computeIfAbsent(collapsedKey, key -> Maps.newTreeMap());
            mFixedLowerGroupIndices.get(collapsedKey).computeIfAbsent(coords.PositionLower, key -> new OneDGridMap<>());
            mFixedLowerGroupIndices.get(collapsedKey).get(coords.PositionLower).put(coords.PositionUpper, groupIndex);

            mFixedUpperGroupIndices.computeIfAbsent(collapsedKey, key -> Maps.newTreeMap());
            mFixedUpperGroupIndices.get(collapsedKey).computeIfAbsent(coords.PositionUpper, key -> new OneDGridMap<>());
            mFixedUpperGroupIndices.get(collapsedKey).get(coords.PositionUpper).put(coords.PositionLower, groupIndex);

            mFixedLowerCoordIndices.computeIfAbsent(collapsedDuplexKey, key -> Maps.newTreeMap());
            mFixedLowerCoordIndices.get(collapsedDuplexKey).computeIfAbsent(coords.PositionLower, key -> new OneDGridMap<>());
            mFixedLowerCoordIndices.get(collapsedDuplexKey).get(coords.PositionLower).put(coords.PositionUpper, groupIndex);

            mFixedUpperCoordIndices.computeIfAbsent(collapsedDuplexKey, key -> Maps.newTreeMap());
            mFixedUpperCoordIndices.get(collapsedDuplexKey).computeIfAbsent(coords.PositionUpper, key -> new OneDGridMap<>());
            mFixedUpperCoordIndices.get(collapsedDuplexKey).get(coords.PositionUpper).put(coords.PositionLower, groupIndex);
        }

        private void collapseDuplexCoordsMerger()
        {
            for(OneDGridMap<FragmentCoords> coordGroup : mCoordsWithUnmappedReadOrMate.values())
            {
                List<List<FragmentCoords>> partitions = coordGroup.partitionValuesByDistance(SINGLE_END_JITTER_COLLAPSE_DISTANCE);
                for(List<FragmentCoords> partition : partitions)
                {
                    for(int i = 1; i < partition.size(); i++)
                        mDuplexCoordsMerger.merge(partition.get(0), partition.get(i));
                }
            }

            List<List<Integer>> mergedCoordGroupsIndices = mFixedLowerCoordIndices.values().stream()
                    .flatMap(x -> x.values().stream())
                    .flatMap(x -> x.partitionValuesByDistance(SINGLE_END_JITTER_COLLAPSE_DISTANCE).stream())
                    .collect(Collectors.toList());

            for(List<Integer> mergedCoordGroupIndices : mergedCoordGroupsIndices)
            {
                FragmentCoords firstCoord = mFullyMappedGroups.get(mergedCoordGroupIndices.get(0)).fragmentCoordinates();
                for(int i = 1; i < mergedCoordGroupIndices.size(); i++)
                {
                    FragmentCoords secondCoord = mFullyMappedGroups.get(mergedCoordGroupIndices.get(i)).fragmentCoordinates();
                    mDuplexCoordsMerger.merge(firstCoord, secondCoord);
                }
            }

            mergedCoordGroupsIndices = mFixedUpperCoordIndices.values().stream()
                    .flatMap(x -> x.values().stream())
                    .flatMap(x -> x.partitionValuesByDistance(SINGLE_END_JITTER_COLLAPSE_DISTANCE).stream())
                    .collect(Collectors.toList());

            for(List<Integer> mergedCoordGroupIndices : mergedCoordGroupsIndices)
            {
                FragmentCoords firstCoord = mFullyMappedGroups.get(mergedCoordGroupIndices.get(0)).fragmentCoordinates();
                for(int i = 1; i < mergedCoordGroupIndices.size(); i++)
                {
                    FragmentCoords secondCoord = mFullyMappedGroups.get(mergedCoordGroupIndices.get(i)).fragmentCoordinates();
                    mDuplexCoordsMerger.merge(firstCoord, secondCoord);
                }
            }
        }

        private List<DuplicateGroup> replaceCoordsWithRepresentatives(final List<DuplicateGroup> collapsedGroups)
        {
            List<DuplicateGroup> groupsWithRepCoord = Lists.newArrayList();
            for(DuplicateGroup collapsedGroup : collapsedGroups)
            {
                FragmentCoords coord = collapsedGroup.fragmentCoordinates();
                FragmentCoords repCoord = mDuplexCoordsMerger.getRepresentative(coord);
                FragmentCoords collapsedCoord = repCoord.withFragmentOrientation(coord.FragmentOrient);

                List<SAMRecord> reads = collapsedGroup.reads();
                groupsWithRepCoord.add(new DuplicateGroup(reads, collapsedCoord));
            }

            return groupsWithRepCoord;
        }

        public FragmentCoordReads getCollapsedGroups()
        {
            if(mFullyMappedGroups.isEmpty() && mGroupsWithUnmappedReadOrMate.isEmpty())
                return null;

            collapseDuplexCoordsMerger();

            List<DuplicateGroup> finalCollapsedGroups = Lists.newArrayList();
            for(OneDGridMap<DuplicateGroup> duplicateGroups : mGroupsWithUnmappedReadOrMate.values())
            {
                finalCollapsedGroups.addAll(
                        duplicateGroups.mergeValuesByDistance(SINGLE_END_JITTER_COLLAPSE_DISTANCE, DUPLICATE_GROUP_MERGER));
            }

            UnionFind<Integer> merger = new UnionFind<>();
            for(int i = 0; i < mFullyMappedGroups.size(); i++)
                merger.add(i);

            List<List<Integer>> mergedGroupsIndices = mFixedLowerGroupIndices.values().stream()
                    .flatMap(x -> x.values().stream())
                    .flatMap(x -> x.partitionValuesByDistance(SINGLE_END_JITTER_COLLAPSE_DISTANCE).stream())
                    .collect(Collectors.toList());

            for(List<Integer> mergedGroupIndices : mergedGroupsIndices)
            {
                for(int i = 1; i < mergedGroupIndices.size(); i++)
                    merger.merge(mergedGroupIndices.get(0), mergedGroupIndices.get(i));
            }

            mergedGroupsIndices = mFixedUpperGroupIndices.values().stream()
                    .flatMap(x -> x.values().stream())
                    .flatMap(x -> x.partitionValuesByDistance(SINGLE_END_JITTER_COLLAPSE_DISTANCE).stream())
                    .collect(Collectors.toList());

            for(List<Integer> mergedGroupIndices : mergedGroupsIndices)
            {
                for(int i = 1; i < mergedGroupIndices.size(); i++)
                    merger.merge(mergedGroupIndices.get(0), mergedGroupIndices.get(i));
            }

            for(Set<Integer> mergedGroupIndices : merger.getPartitions())
            {
                finalCollapsedGroups.add(mergedGroupIndices.stream()
                        .map(mFullyMappedGroups::get)
                        .reduce(DUPLICATE_GROUP_MERGER)
                        .orElse(null));
            }

            return getFragmentCoordReads(replaceCoordsWithRepresentatives(finalCollapsedGroups));
        }
    }

    static FragmentCoordReads illuminaCollapse(
            @Nullable final List<DuplicateGroup> duplicateGroups, @Nullable final List<ReadInfo> singleReads)
    {
        boolean unpairedReads;
        if(singleReads != null && !singleReads.isEmpty())
            unpairedReads = singleReads.get(0).coordinates().Unpaired;
        else
            unpairedReads = duplicateGroups.get(0).fragmentCoordinates().Unpaired;

        if(unpairedReads)
            return new FragmentCoordReads(duplicateGroups, singleReads);

        IlluminaCollapser collapser = new IlluminaCollapser();

        if(singleReads != null)
            singleReads.forEach(collapser::addSingleRead);

        if(duplicateGroups != null)
            duplicateGroups.forEach(collapser::addDuplicateGroup);

        return collapser.getCollapsedGroups();
    }

    class UltimaCollapser
    {
        private final Map<String, OneDGridMap<DuplicateGroup>> mFivePrimeGroups;

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
            OneDGridMap<DuplicateGroup> fivePrimeGroup = getOrCreateFivePrimeGroup(fivePrimeKey);
            fivePrimeGroup.merge(fragEndPos, duplicateGroup, DUPLICATE_GROUP_MERGER);
        }

        public FragmentCoordReads getCollapsedGroups()
        {
            if(mFivePrimeGroups.isEmpty())
                return null;

            List<DuplicateGroup> finalCollapsedGroups = null;
            for(OneDGridMap<DuplicateGroup> fivePrimeGroup : mFivePrimeGroups.values())
            {
                List<DuplicateGroup> collapsedGroups = fivePrimeGroup.mergeValuesByDistance(
                        SINGLE_END_JITTER_COLLAPSE_DISTANCE, DUPLICATE_GROUP_MERGER);

                if(finalCollapsedGroups == null)
                {
                    finalCollapsedGroups = collapsedGroups;
                    continue;
                }

                finalCollapsedGroups.addAll(collapsedGroups);
            }

            return getFragmentCoordReads(finalCollapsedGroups);
        }

        private OneDGridMap<DuplicateGroup> getOrCreateFivePrimeGroup(final String fivePrimeKey)
        {
            mFivePrimeGroups.computeIfAbsent(fivePrimeKey, key -> new OneDGridMap<>());
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
        private final int mMaxDuplicateDistance;
        private final Map<String, TwoDGridMap<DuplicateGroup>> mKeyGroups;

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
            mKeyGroups.computeIfAbsent(collapsedKey, key -> new TwoDGridMap<>());
            mKeyGroups.get(collapsedKey).merge(fragStartPos, fragEndPos, duplicateGroup, DUPLICATE_GROUP_MERGER);
        }

        public FragmentCoordReads getCollapsedGroups()
        {
            if(mKeyGroups.isEmpty())
                return null;

            List<DuplicateGroup> collapsedGroups = Lists.newArrayList();
            for(TwoDGridMap<DuplicateGroup> keyGroup : mKeyGroups.values())
            {
                List<DuplicateGroup> keyGroupCollapsed = keyGroup.mergeValuesByDistance(mMaxDuplicateDistance, DUPLICATE_GROUP_MERGER);

                collapsedGroups.addAll(keyGroupCollapsed);
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
        return new DuplicateGroupCollapser()
        {
            @Override
            public FragmentCoordReads collapse(
                    @Nullable final List<DuplicateGroup> duplicateGroups, @Nullable final List<ReadInfo> singleReads)
            {
                SbxCollapser collapser = new SbxCollapser(maxDuplicateDistance);

                if(singleReads != null)
                    singleReads.forEach(collapser::addSingleRead);

                if(duplicateGroups != null)
                    duplicateGroups.forEach(collapser::addDuplicateGroup);

                return collapser.getCollapsedGroups();
            }
        };
    }
}

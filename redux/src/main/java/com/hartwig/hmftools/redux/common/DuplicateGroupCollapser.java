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
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.collect.OneDGridMap;
import com.hartwig.hmftools.common.collect.TwoDGridMap;
import com.hartwig.hmftools.common.collect.UnionFind;

import org.jetbrains.annotations.Nullable;

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

    class SingleReadOrDuplicateGroup
    {
        private final ReadInfo mSingleRead;
        private final DuplicateGroup mDuplicateGroup;

        public SingleReadOrDuplicateGroup(final ReadInfo singleRead)
        {
            mSingleRead = singleRead;
            mDuplicateGroup = null;
        }

        public SingleReadOrDuplicateGroup(final DuplicateGroup duplicateGroup)
        {
            mSingleRead = null;
            mDuplicateGroup = duplicateGroup;
        }

        public boolean isSingleRead() { return mSingleRead != null; }

        public ReadInfo singleRead() { return mSingleRead; }
        public DuplicateGroup duplicateGroup() { return mDuplicateGroup; }

        public FragmentCoords fragmentCoordinates()
        {
            return mSingleRead != null ? mSingleRead.coordinates() : mDuplicateGroup.fragmentCoordinates();
        }

        public void updateFragmentCoordinates(final FragmentCoords coords)
        {
            if(mSingleRead != null)
                mSingleRead.updateCoordinates(coords);
            else
                mDuplicateGroup.updateFragmentCoordinates(coords);
        }

        public DuplicateGroup asDuplicateGroup()
        {
            if(mDuplicateGroup != null)
                return mDuplicateGroup;

            FragmentCoords preCollapsedCoords = mSingleRead.preCollapsedCoordinates();
            FragmentCoords coords = mSingleRead.coordinates();
            DuplicateGroup duplicateGroup = new DuplicateGroup(Lists.newArrayList(mSingleRead.read()), preCollapsedCoords);
            duplicateGroup.updateFragmentCoordinates(coords);
            return duplicateGroup;
        }

        public SingleReadOrDuplicateGroup merge(final SingleReadOrDuplicateGroup otherGroup)
        {
            return new SingleReadOrDuplicateGroup(asDuplicateGroup().merge(otherGroup.asDuplicateGroup()));
        }
    }

    private static FragmentCoordReads getFragmentCoordReads(final Collection<SingleReadOrDuplicateGroup> collapsedGroups)
    {
        List<DuplicateGroup> duplicateGroups = Lists.newArrayList();
        List<ReadInfo> singleReads = Lists.newArrayList();
        for(SingleReadOrDuplicateGroup collapsedGroup : collapsedGroups)
        {
            if(collapsedGroup.isSingleRead())
            {
                singleReads.add(collapsedGroup.singleRead());
                continue;
            }

            duplicateGroups.add(collapsedGroup.duplicateGroup());
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
    static String collapseToKeyWithoutCoordinates(final FragmentCoords fragmentCoords)
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
        String fragmentOrientationSuffix = fragmentCoords.forwardFragment() ? "" : ":N";
        return format("%s:%s:%s:%s:%s%s%s", fragmentCoords.ChromsomeLower, lowerOrientation, fragmentCoords.ChromsomeUpper, upperOrientation, isLowerString, suppSuffix, fragmentOrientationSuffix);
    }

    int SINGLE_END_JITTER_COLLAPSE_DISTANCE = 10;

    class IlluminaCollapser
    {
        private final List<SingleReadOrDuplicateGroup> mFullyMappedGroups;
        private final Map<String, SortedMap<Integer, OneDGridMap<Integer>>> mFixedLowerGroupIndices;
        private final Map<String, SortedMap<Integer, OneDGridMap<Integer>>> mFixedUpperGroupIndices;
        private final Map<String, OneDGridMap<SingleReadOrDuplicateGroup>> mGroupsWithUnmappedReadOrMate;

        private final Map<String, SortedMap<Integer, OneDGridMap<Integer>>> mFixedLowerForwardCoordIndices;
        private final Map<String, SortedMap<Integer, OneDGridMap<Integer>>> mFixedUpperForwardCoordIndices;
        private final Map<String, OneDGridMap<FragmentCoords>> mForwardCoordsWithUnmappedReadOrMate;

        private final UnionFind<FragmentCoords> mForwardCoordsMerger;

        public IlluminaCollapser()
        {
            mFullyMappedGroups = Lists.newArrayList();
            mFixedLowerGroupIndices = Maps.newHashMap();
            mFixedUpperGroupIndices = Maps.newHashMap();
            mGroupsWithUnmappedReadOrMate = Maps.newHashMap();

            mFixedLowerForwardCoordIndices = Maps.newHashMap();
            mFixedUpperForwardCoordIndices = Maps.newHashMap();
            mForwardCoordsWithUnmappedReadOrMate = Maps.newHashMap();

            mForwardCoordsMerger = new UnionFind<>();
        }

        public void addSingleRead(final ReadInfo readInfo)
        {
            addGroup(new SingleReadOrDuplicateGroup(readInfo));
        }

        public void addDuplicateGroup(final DuplicateGroup duplicateGroup)
        {
            addGroup(new SingleReadOrDuplicateGroup(duplicateGroup));
        }

        private void addGroup(final SingleReadOrDuplicateGroup group)
        {
            FragmentCoords coords = group.fragmentCoordinates();
            FragmentCoords forwardCoords = coords.withFragmentOrientation(FORWARD);
            mForwardCoordsMerger.add(forwardCoords);
            String collapsedKey = collapseToKeyWithoutCoordinates(coords);
            String collapsedForwardKey = collapseToKeyWithoutCoordinates(forwardCoords);
            if(coords.PositionUpper == NO_POSITION)
            {
                mGroupsWithUnmappedReadOrMate.computeIfAbsent(collapsedKey, key -> new OneDGridMap<>());
                mGroupsWithUnmappedReadOrMate.get(collapsedKey).put(coords.PositionLower, group);

                mForwardCoordsWithUnmappedReadOrMate.computeIfAbsent(collapsedForwardKey, key -> new OneDGridMap<>());
                mForwardCoordsWithUnmappedReadOrMate.get(collapsedForwardKey).put(forwardCoords.PositionLower, forwardCoords);
                return;
            }

            mFullyMappedGroups.add(group);
            int groupIndex = mFullyMappedGroups.size() - 1;

            mFixedLowerGroupIndices.computeIfAbsent(collapsedKey, key -> Maps.newTreeMap());
            mFixedLowerGroupIndices.get(collapsedKey).computeIfAbsent(coords.PositionLower, key -> new OneDGridMap<>());
            mFixedLowerGroupIndices.get(collapsedKey).get(coords.PositionLower).put(coords.PositionUpper, groupIndex);

            mFixedUpperGroupIndices.computeIfAbsent(collapsedKey, key -> Maps.newTreeMap());
            mFixedUpperGroupIndices.get(collapsedKey).computeIfAbsent(coords.PositionUpper, key -> new OneDGridMap<>());
            mFixedUpperGroupIndices.get(collapsedKey).get(coords.PositionUpper).put(coords.PositionLower, groupIndex);

            mFixedLowerForwardCoordIndices.computeIfAbsent(collapsedForwardKey, key -> Maps.newTreeMap());
            mFixedLowerForwardCoordIndices.get(collapsedForwardKey).computeIfAbsent(forwardCoords.PositionLower, key -> new OneDGridMap<>());
            mFixedLowerForwardCoordIndices.get(collapsedForwardKey).get(forwardCoords.PositionLower).put(forwardCoords.PositionUpper, groupIndex);

            mFixedUpperForwardCoordIndices.computeIfAbsent(collapsedForwardKey, key -> Maps.newTreeMap());
            mFixedUpperForwardCoordIndices.get(collapsedForwardKey).computeIfAbsent(forwardCoords.PositionUpper, key -> new OneDGridMap<>());
            mFixedUpperForwardCoordIndices.get(collapsedForwardKey).get(forwardCoords.PositionUpper).put(forwardCoords.PositionLower, groupIndex);
        }

        private void collapseForwardCoordsMerger()
        {
            for(OneDGridMap<FragmentCoords> coordGroup : mForwardCoordsWithUnmappedReadOrMate.values())
            {
                List<List<FragmentCoords>> partitions = coordGroup.partitionValuesByDistance(SINGLE_END_JITTER_COLLAPSE_DISTANCE);
                for(List<FragmentCoords> partition : partitions)
                {
                    for(int i = 1; i < partition.size(); i++)
                        mForwardCoordsMerger.merge(partition.get(0), partition.get(i));
                }
            }

            List<List<Integer>> mergedForwardCoordGroupsIndices = mFixedLowerForwardCoordIndices.values().stream()
                    .flatMap(x -> x.values().stream())
                    .flatMap(x -> x.partitionValuesByDistance(SINGLE_END_JITTER_COLLAPSE_DISTANCE).stream())
                    .collect(Collectors.toList());

            for(List<Integer> mergedForwardCoordGroupIndices : mergedForwardCoordGroupsIndices)
            {
                FragmentCoords firstCoord = mFullyMappedGroups.get(mergedForwardCoordGroupIndices.get(0)).fragmentCoordinates();
                FragmentCoords firstForwardCoord = firstCoord.withFragmentOrientation(FORWARD);
                for(int i = 1; i < mergedForwardCoordGroupIndices.size(); i++)
                {
                    FragmentCoords secondCoord = mFullyMappedGroups.get(mergedForwardCoordGroupIndices.get(i)).fragmentCoordinates();
                    FragmentCoords secondForwardCoord = secondCoord.withFragmentOrientation(FORWARD);
                    mForwardCoordsMerger.merge(firstForwardCoord, secondForwardCoord);
                }
            }

            mergedForwardCoordGroupsIndices = mFixedUpperForwardCoordIndices.values().stream()
                    .flatMap(x -> x.values().stream())
                    .flatMap(x -> x.partitionValuesByDistance(SINGLE_END_JITTER_COLLAPSE_DISTANCE).stream())
                    .collect(Collectors.toList());

            for(List<Integer> mergedForwardCoordGroupIndices : mergedForwardCoordGroupsIndices)
            {
                FragmentCoords firstCoord = mFullyMappedGroups.get(mergedForwardCoordGroupIndices.get(0)).fragmentCoordinates();
                FragmentCoords firstForwardCoord = firstCoord.withFragmentOrientation(FORWARD);
                for(int i = 1; i < mergedForwardCoordGroupIndices.size(); i++)
                {
                    FragmentCoords secondCoord = mFullyMappedGroups.get(mergedForwardCoordGroupIndices.get(i)).fragmentCoordinates();
                    FragmentCoords secondForwardCoord = secondCoord.withFragmentOrientation(FORWARD);
                    mForwardCoordsMerger.merge(firstForwardCoord, secondForwardCoord);
                }
            }
        }

        private void replaceCoordsWithRepresentatives(final List<SingleReadOrDuplicateGroup> collapsedGroups)
        {
            for(SingleReadOrDuplicateGroup collapsedGroup : collapsedGroups)
            {
                FragmentCoords coord = collapsedGroup.fragmentCoordinates();
                FragmentCoords forwardCoord = coord.withFragmentOrientation(FORWARD);
                FragmentCoords repForwardCoord = mForwardCoordsMerger.getRepresentative(forwardCoord);
                FragmentCoords collapsedCoord = repForwardCoord.withFragmentOrientation(coord.FragmentOrient);
                collapsedGroup.updateFragmentCoordinates(collapsedCoord);
            }
        }

        public FragmentCoordReads getCollapsedGroups()
        {
            if(mFullyMappedGroups.isEmpty() && mGroupsWithUnmappedReadOrMate.isEmpty())
                return null;

            collapseForwardCoordsMerger();

            List<SingleReadOrDuplicateGroup> finalCollapsedGroups = Lists.newArrayList();
            for(OneDGridMap<SingleReadOrDuplicateGroup> duplicateGroups : mGroupsWithUnmappedReadOrMate.values())
            {
                finalCollapsedGroups.addAll(
                        duplicateGroups.mergeValuesByDistance(SINGLE_END_JITTER_COLLAPSE_DISTANCE, SingleReadOrDuplicateGroup::merge));
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
                        .reduce(SingleReadOrDuplicateGroup::merge)
                        .orElse(null));
            }

            replaceCoordsWithRepresentatives(finalCollapsedGroups);
            return getFragmentCoordReads(finalCollapsedGroups);
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

        var result = collapser.getCollapsedGroups();
        return result;
    }

    class UltimaCollapser
    {
        private final Map<String, OneDGridMap<SingleReadOrDuplicateGroup>> mFivePrimeGroups;

        public UltimaCollapser()
        {
            mFivePrimeGroups = Maps.newHashMap();
        }

        public void addSingleRead(final ReadInfo readInfo)
        {
            addGroup(new SingleReadOrDuplicateGroup(readInfo));
        }

        public void addDuplicateGroup(final DuplicateGroup duplicateGroup)
        {
            addGroup(new SingleReadOrDuplicateGroup(duplicateGroup));
        }

        private void addGroup(final SingleReadOrDuplicateGroup group)
        {
            FragmentCoords coords = group.fragmentCoordinates();
            int fragEndPos = coords.ReadIsLower ? coords.PositionUpper : coords.PositionLower;
            String fivePrimeKey = collapseToFivePrimeKey(coords);
            OneDGridMap<SingleReadOrDuplicateGroup> fivePrimeGroup = getOrCreateFivePrimeGroup(fivePrimeKey);
            fivePrimeGroup.merge(fragEndPos, group, SingleReadOrDuplicateGroup::merge);
        }

        public FragmentCoordReads getCollapsedGroups()
        {
            if(mFivePrimeGroups.isEmpty())
                return null;

            List<SingleReadOrDuplicateGroup> finalCollapsedGroups = null;
            for(OneDGridMap<SingleReadOrDuplicateGroup> fivePrimeGroup : mFivePrimeGroups.values())
            {
                List<SingleReadOrDuplicateGroup> collapsedGroups = fivePrimeGroup.mergeValuesByDistance(
                        SINGLE_END_JITTER_COLLAPSE_DISTANCE, SingleReadOrDuplicateGroup::merge);

                if(finalCollapsedGroups == null)
                {
                    finalCollapsedGroups = collapsedGroups;
                    continue;
                }

                finalCollapsedGroups.addAll(collapsedGroups);
            }

            return getFragmentCoordReads(finalCollapsedGroups);
        }

        private OneDGridMap<SingleReadOrDuplicateGroup> getOrCreateFivePrimeGroup(final String fivePrimeKey)
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
        Map<String, SingleReadOrDuplicateGroup> fivePrimeGroups = Maps.newHashMap();

        if(singleReads != null)
        {
            for(ReadInfo readInfo : singleReads)
            {
                FragmentCoords coords = readInfo.coordinates();
                String fivePrimeKey = collapseToFivePrimeKey(coords);
                fivePrimeGroups.merge(fivePrimeKey, new SingleReadOrDuplicateGroup(readInfo), SingleReadOrDuplicateGroup::merge);
            }
        }

        if(duplicateGroups != null)
        {
            for(DuplicateGroup duplicateGroup : duplicateGroups)
            {
                String fivePrimeKey = collapseToFivePrimeKey(duplicateGroup.fragmentCoordinates());
                fivePrimeGroups.merge(fivePrimeKey, new SingleReadOrDuplicateGroup(duplicateGroup), SingleReadOrDuplicateGroup::merge);
            }
        }

        if(fivePrimeGroups.isEmpty())
            return null;

        return getFragmentCoordReads(fivePrimeGroups.values());
    }

    class SbxCollapser
    {
        private final int mMaxDuplicateDistance;
        private final Map<String, TwoDGridMap<SingleReadOrDuplicateGroup>> mKeyGroups;

        public SbxCollapser(int maxDuplicateDistance)
        {
            mMaxDuplicateDistance = maxDuplicateDistance;
            mKeyGroups = Maps.newHashMap();
        }

        public void addSingleRead(final ReadInfo readInfo)
        {
            addGroup(new SingleReadOrDuplicateGroup(readInfo));
        }

        public void addDuplicateGroup(final DuplicateGroup duplicateGroup)
        {
            addGroup(new SingleReadOrDuplicateGroup(duplicateGroup));
        }

        private void addGroup(final SingleReadOrDuplicateGroup group)
        {
            FragmentCoords coords = group.fragmentCoordinates();
            String collapsedKey = collapseKey(coords);
            int fragStartPos = coords.ReadIsLower ? coords.PositionLower : coords.PositionUpper;
            int fragEndPos = coords.ReadIsLower ? coords.PositionUpper : coords.PositionLower;
            mKeyGroups.computeIfAbsent(collapsedKey, key -> new TwoDGridMap<>());
            mKeyGroups.get(collapsedKey).merge(fragStartPos, fragEndPos, group, SingleReadOrDuplicateGroup::merge);
        }

        public FragmentCoordReads getCollapsedGroups()
        {
            if(mKeyGroups.isEmpty())
                return null;

            List<SingleReadOrDuplicateGroup> collapsedGroups = Lists.newArrayList();
            for(TwoDGridMap<SingleReadOrDuplicateGroup> keyGroup : mKeyGroups.values())
            {
                List<SingleReadOrDuplicateGroup> keyGroupCollapsed = keyGroup.mergeValuesByDistance(
                        mMaxDuplicateDistance, SingleReadOrDuplicateGroup::merge);
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

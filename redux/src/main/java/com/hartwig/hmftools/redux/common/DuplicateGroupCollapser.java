package com.hartwig.hmftools.redux.common;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.sequencing.SequencingType.BIOMODAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ILLUMINA;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ULTIMA;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BinaryOperator;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.collect.OneDGridMap;
import com.hartwig.hmftools.common.collect.TwoDGridMap;
import com.hartwig.hmftools.common.collect.UnionFind;

import org.apache.commons.lang3.NotImplementedException;
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
                ReadInfo singleRead = new ReadInfo(collapsedGroup.reads().get(0), collapsedGroup.fragmentCoordinates_());
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
    static String collapseToKeyWithoutCoordinates(final FragmentCoords fragmentCoords)
    {
        if(fragmentCoords.Unpaired)
            return null;

        if(fragmentCoords.PositionUpper == NO_POSITION)
        {
            String keyWithoutChromosome = fragmentCoords.Key.split(":", 2)[1];
            String[] components = keyWithoutChromosome.split("[:_]", 2);
            if(components.length == 1)
                return "F";

            return components[1].charAt(0) == 'R' ? components[1] : "F_" + components[1];
        }

        String[] components = fragmentCoords.Key.split("_", 3);
        String lowerCoord = components[0];
        String upperCoord = components[1];
        String suffix = components[2];

        String[] lowerCoordComponents = lowerCoord.split(":", 3);
        String lowerOrientation = lowerCoordComponents.length < 3 ? "F" : lowerCoordComponents[2];

        String[] upperCoordComponents = upperCoord.split(":", 3);
        String upperOrientation = upperCoordComponents.length < 3 ? "F" : upperCoordComponents[2];

        return format("%s_%s_%s", lowerOrientation, upperOrientation, suffix);
    }

    int SINGLE_END_JITTER_COLLAPSE_DISTANCE = 10;

    class IlluminaCollapser
    {
        private final Map<String, List<DuplicateGroup>> mCollapedKeyGroups;

        public IlluminaCollapser()
        {
            mCollapedKeyGroups = Maps.newHashMap();
        }

        public void addSingleRead(final ReadInfo readInfo)
        {
            addDuplicateGroup(new DuplicateGroup(null, readInfo.read(), readInfo.coordinates()));
        }

        public void addDuplicateGroup(final DuplicateGroup duplicateGroup)
        {
            FragmentCoords coords = duplicateGroup.fragmentCoordinates_();
            String collpsedKey = collapseToKeyWithoutCoordinates(coords);
            mCollapedKeyGroups.computeIfAbsent(collpsedKey, key -> Lists.newArrayList());
            mCollapedKeyGroups.get(collpsedKey).add(duplicateGroup);
        }

        public FragmentCoordReads getCollapsedGroups()
        {
            if(mCollapedKeyGroups.isEmpty())
                return null;

            List<DuplicateGroup> finalCollapsedGroups = Lists.newArrayList();
            for(List<DuplicateGroup> groups : mCollapedKeyGroups.values())
            {
                UnionFind<Integer> groupMerger = new UnionFind<>();
                for(int i = 0; i < groups.size(); i++)
                    groupMerger.add(i);

                for(int i = 0; i < groups.size() - 1; i++)
                {
                    FragmentCoords coords1 = groups.get(i).fragmentCoordinates_();
                    for(int j = i + 1; j < groups.size(); j++)
                    {
                        FragmentCoords coords2 = groups.get(j).fragmentCoordinates_();
                        if(!coords1.ChromsomeLower.equals(coords2.ChromsomeLower))
                            continue;

                        if(!coords1.ChromsomeUpper.equals(coords2.ChromsomeUpper))
                            continue;

                        int lowerDistance = abs(coords1.PositionLower - coords2.PositionLower);
                        int upperDistance = abs(coords1.PositionUpper - coords2.PositionUpper);
                        if(lowerDistance == 0 && upperDistance <= SINGLE_END_JITTER_COLLAPSE_DISTANCE)
                        {
                            groupMerger.merge(i, j);
                            continue;
                        }

                        if(lowerDistance <= SINGLE_END_JITTER_COLLAPSE_DISTANCE && upperDistance == 0)
                            groupMerger.merge(i, j);
                    }
                }

                Collection<Set<Integer>> mergedGroupsIndices = groupMerger.getPartitions();
                for(Set<Integer> mergedGroupIndices : mergedGroupsIndices)
                {
                    List<DuplicateGroup> duplicateGroups = mergedGroupIndices.stream().map(i -> groups.get(i)).collect(Collectors.toList());
                    DuplicateGroup firstGroup = duplicateGroups.get(0);
                    finalCollapsedGroups.add(firstGroup);
                    for(int i = 1; i < duplicateGroups.size(); i++)
                        firstGroup.addReads(duplicateGroups.get(i).reads());
                }
            }

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
            unpairedReads = duplicateGroups.get(0).fragmentCoordinates_().Unpaired;

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
            FragmentCoords coords = duplicateGroup.fragmentCoordinates_();
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
                String fivePrimeKey = collapseToFivePrimeKey(duplicateGroup.fragmentCoordinates_());
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
            FragmentCoords coords = duplicateGroup.fragmentCoordinates_();
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

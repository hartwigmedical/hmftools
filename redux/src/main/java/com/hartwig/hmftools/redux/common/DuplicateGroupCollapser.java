package com.hartwig.hmftools.redux.common;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.collect.MergeUtils.clusterMerger;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.sequencing.SequencingType.BIOMODAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ULTIMA;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.BiPredicate;
import java.util.function.BinaryOperator;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

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

    static FragmentCoordReads getFragmentCoordReads(final Stream<DuplicateGroup> collapsedGroups)
    {
        final List<DuplicateGroup> duplicateGroups = Lists.newArrayList();
        final List<ReadInfo> singleReads = Lists.newArrayList();
        collapsedGroups.forEach(collapsedGroup ->
        {
            if(collapsedGroup.readCount() == 1)
            {
                singleReads.add(new ReadInfo(collapsedGroup.reads().get(0), collapsedGroup.fragmentCoordinates()));
                return;
            }

            duplicateGroups.add(collapsedGroup);
        });

        return new FragmentCoordReads(duplicateGroups, singleReads);
    }

    static FragmentCoordReads getFragmentCoordReads(final Collection<DuplicateGroup> collapsedGroups)
    {
        return getFragmentCoordReads(collapsedGroups.stream());
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

    @Nullable
    static String collapseToNonOrientedKeyWithoutCoordinates(final FragmentCoords fragmentCoords)
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
        return format("%s:%s:%s:%s:%s%s",
                fragmentCoords.ChromsomeLower, lowerOrientation, fragmentCoords.ChromsomeUpper, upperOrientation, isLowerString, suppSuffix);
    }

    int SINGLE_END_JITTER_COLLAPSE_DISTANCE = 10;
    Comparator<DuplicateGroup> DUPLICATE_GROUP_COMPARATOR = Comparator.comparingInt(DuplicateGroup::readCount).reversed();

    class UltimaCollapser
    {
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
                collapsedGroups.addAll(clusterMerger(
                        fivePrimeGroup,
                        (x, y) -> abs(x - y) <= SINGLE_END_JITTER_COLLAPSE_DISTANCE,
                        DUPLICATE_GROUP_COMPARATOR,
                        DUPLICATE_GROUP_MERGER,
                        null));
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

            BiPredicate<FragStartEnd, FragStartEnd> canMergeFn = (x, y) -> x.distance(y) <= mMaxDuplicateDistance;
            List<DuplicateGroup> collapsedGroups = Lists.newArrayList();
            for(Map<FragStartEnd, DuplicateGroup> keyGroup : mKeyGroups.values())
                collapsedGroups.addAll(clusterMerger(keyGroup, canMergeFn, DUPLICATE_GROUP_COMPARATOR, DUPLICATE_GROUP_MERGER, null));

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

package com.hartwig.hmftools.redux.common;

import static com.hartwig.hmftools.common.sequencing.SequencingType.BIOMODAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ULTIMA;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.collect.OneDGridMap;
import com.hartwig.hmftools.common.collect.TwoDGridMap;

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

    int ULTIMA_MAX_THREE_PRIME_COLLAPSE_DISTANCE = 10;

    class UltimaCollapser
    {
        private final Map<String, OneDGridMap<MultiCoordsDuplicateGroup>> mFivePrimeGroups;

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
            MultiCoordsDuplicateGroup multiCoordGroup = new MultiCoordsDuplicateGroup(duplicateGroup);

            OneDGridMap<MultiCoordsDuplicateGroup> fivePrimeGroup = getOrCreateFivePrimeGroup(fivePrimeKey);
            fivePrimeGroup.merge(fragEndPos, multiCoordGroup, (oldValue, newValue) -> oldValue.merge(newValue));
        }

        public FragmentCoordReads getCollapsedGroups()
        {
            if(mFivePrimeGroups.isEmpty())
                return null;

            List<MultiCoordsDuplicateGroup> finalCollapsedGroups = null;
            for(OneDGridMap<MultiCoordsDuplicateGroup> fivePrimeGroup : mFivePrimeGroups.values())
            {
                List<MultiCoordsDuplicateGroup> collapsedGroups =
                        fivePrimeGroup.mergeValuesByDistance(ULTIMA_MAX_THREE_PRIME_COLLAPSE_DISTANCE, (oldValue, newValue) -> oldValue.merge(newValue));

                if(finalCollapsedGroups == null)
                {
                    finalCollapsedGroups = collapsedGroups;
                    continue;
                }

                finalCollapsedGroups.addAll(collapsedGroups);
            }

            return MultiCoordsFragmentCoordReads.fromCollapsedGroups(finalCollapsedGroups);
        }

        private OneDGridMap<MultiCoordsDuplicateGroup> getOrCreateFivePrimeGroup(final String fivePrimeKey)
        {
            OneDGridMap<MultiCoordsDuplicateGroup> fivePrimeGroup = mFivePrimeGroups.get(fivePrimeKey);
            if(fivePrimeGroup != null)
                return fivePrimeGroup;

            fivePrimeGroup = new OneDGridMap<>();
            mFivePrimeGroups.put(fivePrimeKey, fivePrimeGroup);
            return fivePrimeGroup;
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
        Map<String, MultiCoordsDuplicateGroup> fivePrimeGroups = Maps.newHashMap();

        if(singleReads != null)
        {
            for(ReadInfo readInfo : singleReads)
            {
                SAMRecord read = readInfo.read();
                FragmentCoords coords = readInfo.coordinates();
                String fivePrimeKey = collapseToFivePrimeKey(coords);
                MultiCoordsDuplicateGroup fivePrimeGroup = fivePrimeGroups.get(fivePrimeKey);
                if(fivePrimeGroup == null)
                {
                    fivePrimeGroup = new MultiCoordsDuplicateGroup(read, coords);
                    fivePrimeGroups.put(fivePrimeKey, fivePrimeGroup);
                    continue;
                }

                fivePrimeGroup.addRead(read, coords);
            }
        }

        if(duplicateGroups != null)
        {
            for(DuplicateGroup duplicateGroup : duplicateGroups)
            {
                String fivePrimeKey = collapseToFivePrimeKey(duplicateGroup.fragmentCoordinates());
                MultiCoordsDuplicateGroup fivePrimeGroup = fivePrimeGroups.get(fivePrimeKey);
                if(fivePrimeGroup == null)
                {
                    fivePrimeGroup = new MultiCoordsDuplicateGroup(duplicateGroup);
                    fivePrimeGroups.put(fivePrimeKey, fivePrimeGroup);
                    continue;
                }

                fivePrimeGroup.merge(duplicateGroup);
            }
        }

        if(fivePrimeGroups.isEmpty())
            return null;

        return MultiCoordsFragmentCoordReads.fromCollapsedGroups(fivePrimeGroups.values());
    }

    class SbxCollapser
    {
        private final int mMaxDuplicateDistance;
        private final Map<String, TwoDGridMap<MultiCoordsDuplicateGroup>> mKeyGroups;

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

            TwoDGridMap<MultiCoordsDuplicateGroup> keyGroup = mKeyGroups.get(collapsedKey);
            if(keyGroup == null)
            {
                keyGroup = new TwoDGridMap<>();
                mKeyGroups.put(collapsedKey, keyGroup);
            }

            MultiCoordsDuplicateGroup multiCoordsDuplicateGroup = new MultiCoordsDuplicateGroup(duplicateGroup);
            keyGroup.merge(fragStartPos, fragEndPos, multiCoordsDuplicateGroup, (acc, x) -> acc.merge(x));
        }

        public FragmentCoordReads getCollapsedGroups()
        {
            if(mKeyGroups.isEmpty())
                return null;

            List<MultiCoordsDuplicateGroup> collapsedGroups = Lists.newArrayList();
            for(TwoDGridMap<MultiCoordsDuplicateGroup> keyGroup : mKeyGroups.values())
            {
                List<MultiCoordsDuplicateGroup> keyGroupCollapsed =
                        keyGroup.mergeValuesByDistance(mMaxDuplicateDistance, (acc, x) -> acc.merge(x));

                collapsedGroups.addAll(keyGroupCollapsed);
            }

            return MultiCoordsFragmentCoordReads.fromCollapsedGroups(collapsedGroups);
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

package com.hartwig.hmftools.redux.common;

import static com.hartwig.hmftools.common.sequencing.SequencingType.ULTIMA;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.collect.OneDGridMap;
import com.hartwig.hmftools.common.sequencing.SequencingType;

import org.jetbrains.annotations.Nullable;

@FunctionalInterface
public interface DuplicateGroupCollapser
{
    FragmentCoordReads collapse(@Nullable final List<DuplicateGroup> duplicateGroups, @Nullable final List<ReadInfo> singleReads);

    static DuplicateGroupCollapser fromSequencingType(final SequencingType sequencingType)
    {
        if(sequencingType == ULTIMA)
            return DuplicateGroupCollapser::ultimaCollapse;

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
}

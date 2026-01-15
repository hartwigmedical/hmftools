package com.hartwig.hmftools.redux.duplicate;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.redux.ReduxConstants.SINGLE_END_JITTER_COLLAPSE_DISTANCE;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.redux.common.ReadInfo;

public class JitterUmiMerger
{
    private final UmiConfig mUmiConfig;

    public JitterUmiMerger(final UmiConfig umiConfig)
    {
        mUmiConfig = umiConfig;
    }

    // TODO: form a base class for these merge groups across poly-G, jitter and SBX-collapsing?

    private class JitterMergeGroup
    {
        public boolean KeyOnLower;
        public final FragmentCoords FragmentCoords;
        public final List<Object> Fragments;
        public boolean Merged;

        public JitterMergeGroup(final Object fragment, final FragmentCoords fragmentCoords, final boolean keyOnLower)
        {
            KeyOnLower = keyOnLower;
            FragmentCoords = fragmentCoords;
            Fragments = Lists.newArrayList(fragment);
            Merged = false;
        }

        public boolean canMerge(final JitterMergeGroup other)
        {
            if(Merged || other.Merged)
                return false;

            if(KeyOnLower)
            {
                // other coords must have same orientation and match within the permitted range
                return FragmentCoords.OrientUpper == other.FragmentCoords.OrientUpper
                    && abs(FragmentCoords.PositionUpper - other.FragmentCoords.PositionLower) <= SINGLE_END_JITTER_COLLAPSE_DISTANCE;
            }
            else
            {
                return FragmentCoords.OrientLower == other.FragmentCoords.OrientLower
                    && abs(FragmentCoords.PositionLower - other.FragmentCoords.PositionLower) <= SINGLE_END_JITTER_COLLAPSE_DISTANCE;
            }
        }

        public void merge(final JitterMergeGroup other)
        {
            Fragments.addAll(other.Fragments);
            other.Merged = true;
        }

        public String toString()
        {
            return format("coords(%s) frags(%d) merged(%s)", FragmentCoords, Fragments.size(), Merged);
        }
    }

    public void mergeGroups(final List<DuplicateGroup> umiGroups, final List<ReadInfo> singleFragments)
    {
        if(umiGroups.isEmpty() && singleFragments.isEmpty())
            return;

        Map<String,List<JitterMergeGroup>> lowerPosGroups = Maps.newHashMap();
        Map<String,List<JitterMergeGroup>> upperPosGroups = Maps.newHashMap();

        // collapse fragments with matching UMIs after trimming and the same 5' unclipped position
        boolean hasCandidates = false;
        for(DuplicateGroup duplicateGroup : umiGroups)
        {
            for(int g = 0; g <= 1; ++g)
            {
                boolean useLower = (g == 0);

                String key = formKey(duplicateGroup.umi(), duplicateGroup.fragmentCoordinates(), useLower);
                Map<String,List<JitterMergeGroup>> posGroups = useLower ? lowerPosGroups : upperPosGroups;

                List<JitterMergeGroup> groups = posGroups.get(key);

                if(groups == null)
                {
                    groups = Lists.newArrayList();
                    posGroups.put(key, groups);
                }
                else
                {
                    hasCandidates = true;
                }

                groups.add(new JitterMergeGroup(duplicateGroup, duplicateGroup.fragmentCoordinates(), useLower));
            }
        }

        for(ReadInfo readInfo : singleFragments)
        {
            String umiId = mUmiConfig.extractUmiId(readInfo.read().getReadName());

            for(int g = 0; g <= 1; ++g)
            {
                boolean useLower = (g == 0);

                String key = formKey(umiId, readInfo.coordinates(), useLower);
                Map<String,List<JitterMergeGroup>> posGroups = useLower ? lowerPosGroups : upperPosGroups;

                List<JitterMergeGroup> groups = posGroups.get(key);

                if(groups == null)
                {
                    groups = Lists.newArrayList();
                    posGroups.put(key, groups);
                }
                else
                {
                    hasCandidates = true;
                }

                groups.add(new JitterMergeGroup(readInfo, readInfo.coordinates(), useLower));
            }
        }

        if(!hasCandidates)
            return;

        boolean hasMerges = false;

        for(int g = 0; g <= 1; ++g)
        {
            boolean useLower = (g == 0);

            Map<String,List<JitterMergeGroup>> posGroups = useLower ? lowerPosGroups : upperPosGroups;

            for(List<JitterMergeGroup> groups : posGroups.values())
            {
                if(groups.size() == 1)
                    continue;

                for(int i = 0; i < groups.size() - 1; ++i)
                {
                    JitterMergeGroup group1 = groups.get(i);

                    for(int j = i + 1; j < groups.size(); ++j)
                    {
                        JitterMergeGroup group2 = groups.get(j);

                        if(group1.FragmentCoords.PositionLower != group2.FragmentCoords.PositionLower)
                            break;

                        if(group1.canMerge(group2))
                        {
                            group1.merge(group2);
                            hasMerges = true;
                        }
                    }
                }
            }
        }

        if(!hasMerges)
            return;

        for(int g = 0; g <= 1; ++g)
        {
            boolean useLower = (g == 0);

            Map<String,List<JitterMergeGroup>> posGroups = useLower ? lowerPosGroups : upperPosGroups;

            for(List<JitterMergeGroup> groups : posGroups.values())
            {
                if(groups.size() == 1)
                    continue;

                for(JitterMergeGroup group : groups)
                {
                    if(group.Merged)
                    {
                        for(Object fragment : group.Fragments)
                        {
                            if(fragment instanceof DuplicateGroup)
                            {
                                umiGroups.remove(fragment);
                            }
                            else
                            {
                                singleFragments.remove(fragment);
                            }
                        }
                    }
                    else if(group.Fragments.size() > 1)
                    {
                        DuplicateGroup originalGroup = null;

                        for(int i = 0; i < group.Fragments.size(); ++i)
                        {
                            Object fragment = group.Fragments.get(i);

                            ReadInfo readInfo = fragment instanceof ReadInfo ? (ReadInfo) fragment : null;
                            DuplicateGroup duplicateGroup = fragment instanceof DuplicateGroup ? (DuplicateGroup) fragment : null;

                            if(originalGroup == null)
                            {
                                if(readInfo != null)
                                {
                                    // convert to duplicate group
                                    String umiId = mUmiConfig.extractUmiId(readInfo.read().getReadName());
                                    originalGroup = new DuplicateGroup(umiId, readInfo.read(), readInfo.coordinates());
                                    umiGroups.add(originalGroup);
                                    singleFragments.remove(readInfo);
                                }
                                else
                                {
                                    originalGroup = duplicateGroup;
                                }
                            }
                            else
                            {
                                if(readInfo != null)
                                {
                                    originalGroup.addReads(List.of(readInfo.read()));
                                }
                                else
                                {
                                    originalGroup.addReads(duplicateGroup.allReads());
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    private static String formKey(final String umi, final FragmentCoords coords, final boolean keyOnLower)
    {
        // same UMI
        Orientation orient = keyOnLower ? coords.OrientLower : coords.OrientUpper;
        int position = keyOnLower ? coords.PositionLower : coords.PositionUpper;
        return format("%s_%d_%d_%s", umi, position, orient, coords.SuppReadInfo != null);
    }
}

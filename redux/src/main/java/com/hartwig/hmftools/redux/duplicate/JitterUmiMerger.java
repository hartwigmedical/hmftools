package com.hartwig.hmftools.redux.duplicate;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.redux.ReduxConstants.SINGLE_END_JITTER_COLLAPSE_DISTANCE;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_ORIENT_FORWARD;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_ORIENT_REVERSE;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_READ_SUPP_INFO;
import static com.hartwig.hmftools.redux.duplicate.UmiGroupBuilder.hasDuplexUmiMatch;
import static com.hartwig.hmftools.redux.duplicate.UmiUtils.exceedsUmiIdDiff;

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

    private class JitterMergeGroup
    {
        public final String Umi;
        public final FragmentCoords FragmentCoords;
        public final List<Object> Fragments;
        public boolean Merged;

        public JitterMergeGroup(final Object fragment, final FragmentCoords fragmentCoords, final String umi)
        {
            Umi = umi;
            FragmentCoords = fragmentCoords;
            Fragments = Lists.newArrayList(fragment);
            Merged = false;
        }

        public boolean canMerge(final JitterMergeGroup other, boolean checkLower)
        {
            if(other.Merged)
                return false;

            if(FragmentCoords.ReadIsLower != other.FragmentCoords.ReadIsLower)
                return false; // this could be the same fragment (ie same reads) but the lower (R1) trying to merge with the upper (R2) reads

            // check UMIs - exact matches are required where jitter was used
            if(FragmentCoords.FragmentOrient == other.FragmentCoords.FragmentOrient)
            {
                if(exceedsUmiIdDiff(Umi, other.Umi, 0))
                //if(!Umi.equals(other.Umi))
                    return false;
            }
            else
            {
                if(!hasDuplexUmiMatch(Umi, other.Umi, mUmiConfig.DuplexDelim, 0))
                    return false;
            }

            if(checkLower)
            {
                return FragmentCoords.OrientLower == other.FragmentCoords.OrientLower
                    && abs(FragmentCoords.PositionLower - other.FragmentCoords.PositionLower) <= SINGLE_END_JITTER_COLLAPSE_DISTANCE;
            }
            else
            {
                // other coords must have same orientation and match within the permitted range
                return FragmentCoords.OrientUpper == other.FragmentCoords.OrientUpper
                    && abs(FragmentCoords.PositionUpper - other.FragmentCoords.PositionUpper) <= SINGLE_END_JITTER_COLLAPSE_DISTANCE;
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
        // convert each duplicate group and single fragment to a candidate jitter group, and add each of these to 2 separate maps
        // so as to test the lower and upper coordinates in turn being jitter-matched

        boolean hasCandidates = false;
        for(DuplicateGroup duplicateGroup : umiGroups)
        {
            JitterMergeGroup group = new JitterMergeGroup(duplicateGroup, duplicateGroup.fragCoordinates(), duplicateGroup.umi());

            for(int g = 0; g <= 1; ++g)
            {
                boolean useLower = (g == 0);

                String key = formKey(duplicateGroup.fragCoordinates(), useLower);
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

                groups.add(group);
            }
        }

        for(ReadInfo readInfo : singleFragments)
        {
            String umiId = readInfo.getOrExtractUmi(mUmiConfig);
            JitterMergeGroup group = new JitterMergeGroup(readInfo, readInfo.fragCoordinates(), umiId);

            for(int g = 0; g <= 1; ++g)
            {
                boolean useLower = (g == 0);

                String key = formKey(readInfo.fragCoordinates(), useLower);
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

                groups.add(group);
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

                    if(group1.Merged)
                        continue;

                    for(int j = i + 1; j < groups.size(); ++j)
                    {
                        JitterMergeGroup group2 = groups.get(j);

                        if(group1.canMerge(group2, !useLower))
                        {
                            hasMerges = true;

                            if(prioritiseFirstCoords(group1.FragmentCoords, group2.FragmentCoords))
                            {
                                group1.merge(group2);
                            }
                            else
                            {
                                group2.merge(group1);
                                break;
                            }
                        }
                    }
                }
            }
        }

        if(!hasMerges)
            return;

        for(List<JitterMergeGroup> groups : lowerPosGroups.values()) // both maps contains the same groups, just keyed differently
        {
            for(JitterMergeGroup group : groups)
            {
                if(group.Merged)
                {
                    // removed groups single fragments which have been merged
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
                    // and then add in the reads from merged groups and single fragments to their new parent group
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
                                String umiId = readInfo.getOrExtractUmi(mUmiConfig);
                                originalGroup = new DuplicateGroup(umiId, readInfo.read(), readInfo.fragCoordinates());
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
                                originalGroup.addNonConsensusReads(List.of(readInfo.read()));

                                if(originalGroup.fragCoordinates().FragmentOrient != readInfo.fragCoordinates().FragmentOrient)
                                    originalGroup.registerDualStrand();
                            }
                            else
                            {
                                originalGroup.addNonConsensusReads(duplicateGroup.allReads());

                                if(originalGroup.fragCoordinates().FragmentOrient != duplicateGroup.fragCoordinates().FragmentOrient)
                                    originalGroup.registerDualStrand();
                            }
                        }
                    }
                }
            }
        }
    }

    private static String formKey(final FragmentCoords coords, final boolean keyOnLower)
    {
        // same UMI
        Orientation orient = keyOnLower ? coords.OrientLower : coords.OrientUpper;
        int position = keyOnLower ? coords.PositionLower : coords.PositionUpper;

        return format("%d_%c_%c",
                position, orient.isForward() ? COORD_ORIENT_FORWARD : COORD_ORIENT_REVERSE,
                coords.SuppReadInfo != null ? COORD_READ_SUPP_INFO : 'P');
    }

    private static boolean prioritiseFirstCoords(final FragmentCoords coords1, final FragmentCoords coords2)
    {
        return (coords1.PositionLower <= coords2.PositionLower) && (coords1.PositionUpper <= coords2.PositionUpper);
    }
}

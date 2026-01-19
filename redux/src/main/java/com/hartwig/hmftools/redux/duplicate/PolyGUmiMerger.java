package com.hartwig.hmftools.redux.duplicate;

import static java.lang.Math.min;
import static java.lang.Math.sin;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.redux.ReduxConstants.MIN_POLYG_UMI_TAIL_LENGTH;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_ORIENT_FORWARD;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_ORIENT_REVERSE;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_READ_SUPP_INFO;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_READ_UNMAPPED;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.redux.common.ReadInfo;

import htsjdk.samtools.SAMRecord;

public class PolyGUmiMerger
{
    // applies to Illumina with duplex UMIs
    // looks for fragments where the only difference in UMI can be explained by degradation to Gs on the unmapped read side
    // fragments with both reads mapped take precedence for consensus, and in this case the unmapped reads are all marked as duplicates

    private final UmiConfig mUmiConfig;

    public PolyGUmiMerger(final UmiConfig umiConfig)
    {
        mUmiConfig = umiConfig;
    }

    private class PolyGMergeGroup
    {
        public final List<String> Umis;
        public final boolean UmiHasPolyGTail;
        public final FragmentCoords FragmentCoords;
        public final List<Object> Fragments;
        public boolean Merged;

        public PolyGMergeGroup(final String umi, final Object fragment, final FragmentCoords fragmentCoords)
        {
            Umis = Lists.newArrayList(trimPolyGTail(umi));
            UmiHasPolyGTail = hasPolyGTail(umi);

            FragmentCoords = fragmentCoords;
            Fragments = Lists.newArrayList(fragment);
            Merged = false;
        }

        public PolyGMergeGroup(final DuplicateGroup duplicateGroup)
        {
            this(duplicateGroup.umi(), duplicateGroup, duplicateGroup.fragCoordinates());
        }

        public PolyGMergeGroup(final ReadInfo singleRead, final String umi)
        {
            this(umi, singleRead, singleRead.fragCoordinates());
        }

        public boolean canMerge(final PolyGMergeGroup other)
        {
            if(other.Merged)
                return false;

            if(!fragmentCoordinatesMatch(FragmentCoords, other.FragmentCoords))
                return false;

            for(String umi : Umis)
            {
                if(other.Umis.stream().anyMatch(x -> umiMatch(umi, x)))
                    return true;
            }

            return false;
        }

        public void merge(final PolyGMergeGroup other)
        {
            Umis.addAll(other.Umis);
            Fragments.addAll(other.Fragments);
            other.Merged = true;
        }

        public String toString()
        {
            return format("umi(%s) coords(%s) frags(%d) merged(%s)", Umis.get(0), FragmentCoords, Fragments.size(), Merged);
        }
    }

    public void mergeGroups(final List<DuplicateGroup> umiGroups, final List<ReadInfo> singleFragments)
    {
        if(umiGroups.isEmpty() && singleFragments.isEmpty())
            return;

        // first a quick search for either poly-G tails or unmapped reads
        if(!hasCandidateFragments(umiGroups, singleFragments))
            return;

        // at least one fragment must have unmapped mates and a poly-G tail on the UMI
        Map<String,List<PolyGMergeGroup>> groupMap = findCandidateFragments(umiGroups, singleFragments);

        if(groupMap.isEmpty())
            return;

        boolean hasMerges = false;

        for(List<PolyGMergeGroup> groups : groupMap.values())
        {
            if(groups.size() == 1)
                continue;

            for(int i = 0; i < groups.size() - 1; ++i)
            {
                PolyGMergeGroup group1 = groups.get(i);

                if(group1.Merged)
                    continue;

                for(int j = i + 1; j < groups.size(); ++j)
                {
                    PolyGMergeGroup group2 = groups.get(j);

                    if(group1.canMerge(group2))
                    {
                        group1.merge(group2);
                        hasMerges = true;
                    }
                }
            }
        }

        if(!hasMerges)
            return;

        // find groups with mapped remote reads, so as to mark matching unmapped reads
        for(List<PolyGMergeGroup> groups : groupMap.values())
        {
            if(groups.size() == 1)
                continue;

            for(PolyGMergeGroup group : groups)
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
                                originalGroup = new DuplicateGroup(readInfo.umi(), readInfo.read(), readInfo.fragCoordinates());
                                umiGroups.add(originalGroup);
                                group.Fragments.set(0, originalGroup);
                                singleFragments.remove(readInfo);
                            }
                            else
                            {
                                originalGroup = duplicateGroup;
                            }
                        }
                        else
                        {
                            List<SAMRecord> readsToAdd = readInfo != null ? List.of(readInfo.read()) : duplicateGroup.reads();

                            if(originalGroup.fragCoordinates().UnmappedSourced)
                            {
                                originalGroup.addNonConsensusReads(readsToAdd);
                            }
                            else
                            {
                                originalGroup.addReads(readsToAdd);
                            }
                        }
                    }
                }
            }
        }

        // finally mark any reads or groups which contain only unmapped reads
        List<UnmappingKey> groupsWithMappedMates = Lists.newArrayList();

        for(List<PolyGMergeGroup> groups : groupMap.values())
        {
            for(PolyGMergeGroup group : groups)
            {
                if(group.Merged || group.FragmentCoords.UnmappedSourced)
                    continue;

                if(group.Fragments.get(0) instanceof DuplicateGroup)
                {
                    DuplicateGroup duplicateGroup = (DuplicateGroup)group.Fragments.get(0);

                    if(duplicateGroup.reads().stream().anyMatch(x -> !x.getMateUnmappedFlag()))
                    {
                        // use trimmed UMI
                        groupsWithMappedMates.add(new UnmappingKey(group.Umis.get(0), duplicateGroup.fragCoordinates()));
                    }
                }
            }
        }

        for(List<PolyGMergeGroup> groups : groupMap.values())
        {
            for(PolyGMergeGroup group : groups)
            {
                if(group.Merged || !group.FragmentCoords.UnmappedSourced)
                    continue;

                if(groupsWithMappedMates.stream().anyMatch(x -> x.matches(group.Umis.get(0), group.FragmentCoords)))
                {
                    if(group.Fragments.get(0) instanceof DuplicateGroup)
                    {
                        DuplicateGroup duplicateGroup = (DuplicateGroup)group.Fragments.get(0);

                        if(duplicateGroup.reads().stream().noneMatch(x -> x.getMateUnmappedFlag()))
                            continue;

                        duplicateGroup.markPolyGUnmapped();
                    }
                    else
                    {
                        ReadInfo readInfo = (ReadInfo)group.Fragments.get(0);

                        if(!readInfo.read().getReadUnmappedFlag())
                            continue;

                        DuplicateGroup duplicateGroup = new DuplicateGroup(readInfo.umi(), readInfo.read(), readInfo.fragCoordinates());
                        duplicateGroup.markPolyGUnmapped();
                        umiGroups.add(duplicateGroup);
                        singleFragments.remove(readInfo);
                    }
                }
            }
        }
    }

    private static boolean umiMatch(final String umi1, final String umi2)
    {
        for(int i = 0; i < min(umi1.length(), umi2.length()); ++i)
        {
            if(umi1.charAt(i) != umi2.charAt(i))
                return false;
        }

        return true;
    }

    private static String formKey(final FragmentCoords coords)
    {
        return format("%d_%c_%c_%s",
                coords.PositionLower, coords.OrientLower.isForward() ? COORD_ORIENT_FORWARD : COORD_ORIENT_REVERSE,
                coords.UnmappedSourced ? COORD_READ_UNMAPPED : 'M',
                coords.SuppReadInfo != null ? COORD_READ_SUPP_INFO : 'P');
    }

    private static boolean fragmentCoordinatesMatch(final FragmentCoords coords1, final FragmentCoords coords2)
    {
        // require a match on the lower coordinates
        return coords1.PositionLower == coords2.PositionLower
            && coords1.UnmappedSourced == coords2.UnmappedSourced
            && coords1.OrientLower == coords2.OrientLower
            && (coords1.SuppReadInfo != null) == (coords2.SuppReadInfo != null);
    }

    private boolean hasCandidateFragments(final List<DuplicateGroup> umiGroups, final List<ReadInfo> singleFragments)
    {
        boolean hasPolyGUmi = false;
        boolean hasUnmappedReads = false;

        for(DuplicateGroup duplicateGroup : umiGroups)
        {
            hasPolyGUmi |= hasUnmappedRead(duplicateGroup.fragCoordinates());
            hasUnmappedReads |= hasPolyGTail(duplicateGroup.umi());

            if(hasPolyGUmi && hasUnmappedReads)
                break;
        }

        for(ReadInfo readInfo : singleFragments)
        {
            hasPolyGUmi |= hasUnmappedRead(readInfo.fragCoordinates());
            hasUnmappedReads |= hasPolyGTail(readInfo.getOrExtractUmi(mUmiConfig));

            if(hasPolyGUmi && hasUnmappedReads)
                break;
        }

        return hasPolyGUmi && hasUnmappedReads;
    }

    private Map<String,List<PolyGMergeGroup>> findCandidateFragments(
            final List<DuplicateGroup> umiGroups, final List<ReadInfo> singleFragments)
    {
        Map<String,List<PolyGMergeGroup>> groupMap = Maps.newHashMap();

        // collapse fragments with matching UMIs after trimming and the same 5' unclipped position
        for(DuplicateGroup duplicateGroup : umiGroups)
        {
            String key = formKey(duplicateGroup.fragCoordinates());
            List<PolyGMergeGroup> groups = groupMap.get(key);

            if(groups == null)
            {
                groups = Lists.newArrayList();
                groupMap.put(key, groups);
            }

            groups.add(new PolyGMergeGroup(duplicateGroup));
        }

        for(ReadInfo readInfo : singleFragments)
        {
            String umi = readInfo.getOrExtractUmi(mUmiConfig);

            String key = formKey(readInfo.fragCoordinates());
            List<PolyGMergeGroup> groups = groupMap.get(key);

            if(groups == null)
            {
                groups = Lists.newArrayList();
                groupMap.put(key, groups);
            }

            groups.add(new PolyGMergeGroup(readInfo, umi));
        }

        List<String> invalidGroups = Lists.newArrayList();

        for(Map.Entry<String,List<PolyGMergeGroup>> entry : groupMap.entrySet())
        {
            List<PolyGMergeGroup> groups = entry.getValue();

            boolean hasPolyGUmi = groups.stream().anyMatch(x -> x.UmiHasPolyGTail);
            boolean hasUnmappedRead = groups.stream().anyMatch(x -> hasUnmappedRead(x.FragmentCoords));

            if(!hasPolyGUmi && !hasUnmappedRead)
                invalidGroups.add(entry.getKey());
        }

        invalidGroups.forEach(x -> groupMap.remove(x));

        return groupMap;
    }

    private static String trimPolyGTail(final String umiId)
    {
        int tailLength;
        for(tailLength = 0; tailLength < umiId.length(); tailLength++)
        {
            if(umiId.charAt(umiId.length() - 1 - tailLength) != 'G')
                break;
        }

        return umiId.substring(0, umiId.length() - tailLength);
    }

    private boolean hasPolyGTail(final String umiId)
    {
        int tailLength;
        for(tailLength = 0; tailLength < umiId.length(); tailLength++)
        {
            if(umiId.charAt(umiId.length() - 1 - tailLength) != 'G')
                break;
        }

        return tailLength >= MIN_POLYG_UMI_TAIL_LENGTH;
    }

    private static boolean hasUnmappedRead(final FragmentCoords fragmentCoords)
    {
        return !fragmentCoords.Unpaired && fragmentCoords.SuppReadInfo == null && fragmentCoords.PositionUpper == NO_POSITION;
    }

    private class UnmappingKey
    {
        public final FragmentCoords Coords;
        public final String TrimmedUmi;

        public UnmappingKey(final String trimmedUmi, final FragmentCoords coords)
        {
            Coords = coords;
            TrimmedUmi = trimmedUmi;
        }

        private boolean matches(final String umi, final FragmentCoords coords)
        {
            return coords.PositionLower == Coords.PositionLower && coords.OrientLower == Coords.OrientLower && umiMatch(TrimmedUmi, umi);
        }
    }
}

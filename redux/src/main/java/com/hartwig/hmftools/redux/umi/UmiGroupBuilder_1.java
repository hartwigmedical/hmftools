package com.hartwig.hmftools.redux.umi;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.redux.common.Constants.MAX_IMBALANCED_UMI_BASE_DIFF;
import static com.hartwig.hmftools.redux.common.Constants.MAX_IMBALANCED_UMI_COUNT;
import static com.hartwig.hmftools.redux.common.Constants.MIN_POLYG_UMI_TAIL_LENGTH;
import static com.hartwig.hmftools.redux.umi.UmiUtils.exceedsUmiIdDiff;
import static com.hartwig.hmftools.redux.umi.UmiUtils.polyGTailLength;
import static com.hartwig.hmftools.redux.umi.UmiUtils.trimPolyGTail;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.collect.UnionFind;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.redux.common.DuplicateGroup;
import com.hartwig.hmftools.redux.common.FragmentCoords;
import com.hartwig.hmftools.redux.common.ReadInfo;

import htsjdk.samtools.SAMRecord;

public class UmiGroupBuilder_1
{
    private final SequencingType mSequencing;
    private final UmiConfig mUmiConfig;

    public UmiGroupBuilder_1(final SequencingType sequencing, final UmiConfig config)
    {
        mSequencing = sequencing;
        mUmiConfig = config;
    }

    public List<DuplicateGroup> processUmiGroups(
            final List<DuplicateGroup> duplicateGroups, final List<ReadInfo> singleFragments, boolean captureStats)
    {
        // the list of input duplicate group and single reads will typically contain unrelated fragments, but any which are related - ie
        // same coordinates but with different orientations need to be handled together

        // IDEA: if R1 and R2 are present in the same collection of fragments, ideally only establish UMI status once and propagate this
        // to the other set of reads

        // organise groups by their UMIs, applying base-difference collapsing rules
        // UMI stats require evaluation of uncollapsed UMI groups with the same coordinates
        // at the same time organise UMI groups by the coordinates

        boolean formCoordGroups = mUmiConfig.BaseStats || (duplicateGroups.size() + singleFragments.size() > 1);

        List<DuplicateGroup> allUmiGroups = Lists.newArrayList();

        for(DuplicateGroup duplicateGroup : duplicateGroups)
        {
            List<DuplicateGroup> umiGroups = buildUmiGroups(duplicateGroup.fragmentCoordinates(), duplicateGroup.reads(), mUmiConfig);
            allUmiGroups.addAll(umiGroups);
        }

        // add in single fragments
        for(ReadInfo readInfo : singleFragments)
        {
            SAMRecord read = readInfo.read();
            FragmentCoords currentCoords = readInfo.coordinates();
            FragmentCoords preCollapsedCoords = readInfo.preCollapsedCoordinates();
            String umiId = mUmiConfig.extractUmiId(read.getReadName());
            DuplicateGroup duplicateGroup = new DuplicateGroup(umiId, Lists.newArrayList(read), preCollapsedCoords);
            duplicateGroup.updateFragmentCoordinates(currentCoords);
            allUmiGroups.add(duplicateGroup);
        }

        singleFragments.clear();

        // TODO: Uncomment.
//        collapsePolyGDuplexUmis(mSequencing, mUmiConfig, allUmiGroups);

        if(formCoordGroups)
        {
            // add in order of descending by fragment count for non-duplex collapsing
            Collections.sort(allUmiGroups, new UmiUtils.SizeComparator());

            // collapse duplex and single UMIs with opposite orientations
            collapseCoordinateGroup(allUmiGroups);
        }

        List<DuplicateGroup> finalUmiGroups = Lists.newArrayList();

        for(DuplicateGroup umiGroup : allUmiGroups)
        {
            if(umiGroup.readCount() == 1)
            {
                // drop any single fragments
                singleFragments.add(new ReadInfo(umiGroup.reads().get(0), umiGroup.fragmentCoordinates()));
                continue;
            }

            finalUmiGroups.add(umiGroup);
        }

        return finalUmiGroups;
    }

    public static List<DuplicateGroup> buildUmiGroups(final FragmentCoords fragCoords, final List<SAMRecord> reads, final UmiConfig config)
    {
        Map<String,DuplicateGroup> groups = Maps.newHashMap();
        boolean checkDefinedUmis = config.hasDefinedUmis();
        boolean useDefinedUmis = checkDefinedUmis;

        for(SAMRecord read : reads)
        {
            String umiId = config.extractUmiId(read.getReadName());

            if(checkDefinedUmis)
            {
                String definedUmiId = config.matchDefinedUmiId(umiId);
                if(definedUmiId == null)
                {
                    useDefinedUmis = false;
                    checkDefinedUmis = false;
                }
                else
                {
                    umiId = definedUmiId;
                }
            }

            DuplicateGroup group = groups.get(umiId);

            if(group == null)
            {
                groups.put(umiId, new DuplicateGroup(umiId, read, fragCoords));
            }
            else
            {
                group.addRead(read);
            }
        }

        if(useDefinedUmis)
        {
            return groups.values().stream().collect(Collectors.toList());
        }

        // order groups by descending number of fragments
        List<DuplicateGroup> orderedGroups = groups.values().stream().sorted(new UmiUtils.SizeComparator()).collect(Collectors.toList());

        // then apply the directional model, where smaller groups are merged into larger ones
        int i = 0;
        while(i < orderedGroups.size() - 1)
        {
            DuplicateGroup first = orderedGroups.get(i);

            List<DuplicateGroup> cluster = Lists.newArrayList(first);

            int j = i + 1;
            while(j < orderedGroups.size())
            {
                DuplicateGroup second = orderedGroups.get(j);

                boolean merged = false;

                for(DuplicateGroup existing : cluster)
                {
                    if(existing.readCount() >= second.readCount() && !exceedsUmiIdDiff(existing.umiId(), second.umiId(), config.PermittedBaseDiff))
                    {
                        merged = true;
                        break;
                    }
                }

                if(!merged)
                {
                    ++j;
                }
                else
                {
                    orderedGroups.remove(j);
                    cluster.add(second);

                    // restart the search since a newly added group may be close enough to a skipped one
                    j = i + 1;
                }
            }

            for(j = 1; j < cluster.size(); ++j)
            {
                first.addReads(cluster.get(j).reads());
            }

            ++i;
        }

        // run a check allowing collapsing of UMIs with 2-base differences
        if(orderedGroups.size() > 1)
        {
            i = 0;
            while(i < orderedGroups.size())
            {
                DuplicateGroup first = orderedGroups.get(i);

                int j = i + 1;
                while(j < orderedGroups.size())
                {
                    DuplicateGroup second = orderedGroups.get(j);

                    if(!exceedsUmiIdDiff(first.umiId(), second.umiId(), config.PermittedBaseDiff + 1))
                    {
                        first.addReads(second.reads());
                        orderedGroups.remove(j);
                    }
                    else
                    {
                        ++j;
                    }
                }

                ++i;
            }
        }

        // run a check allowing collapsing of UMIs with 4-base differences where significant imbalance exists
        boolean hasLargeGroups = orderedGroups.stream().anyMatch(x -> x.readCount() >= MAX_IMBALANCED_UMI_COUNT);

        if(orderedGroups.size() > 1 && hasLargeGroups)
        {
            i = 0;
            while(i < orderedGroups.size())
            {
                DuplicateGroup first = orderedGroups.get(i);

                int j = i + 1;
                while(j < orderedGroups.size())
                {
                    DuplicateGroup second = orderedGroups.get(j);

                    double maxCountRatio = first.readCount() >= second.readCount() ?
                            first.readCount() / (double)second.readCount() : second.readCount() / (double)first.readCount();

                    if(maxCountRatio >= MAX_IMBALANCED_UMI_COUNT && !exceedsUmiIdDiff(first.umiId(), second.umiId(), MAX_IMBALANCED_UMI_BASE_DIFF))
                    {
                        first.addReads(second.reads());
                        orderedGroups.remove(j);
                    }
                    else
                    {
                        ++j;
                    }
                }

                ++i;
            }
        }

        return orderedGroups;
    }

    private void collapseCoordinateGroup(final List<DuplicateGroup> allUmiGroups)
    {
        if(!mUmiConfig.Duplex)
            return;

        List<DuplicateGroup> finalUmiGroups = allUmiGroups.stream().filter(x -> x.fragmentCoordinates().PositionUpper == NO_POSITION).collect(Collectors.toList());
        List<DuplicateGroup> filteredUmiGroups = allUmiGroups.stream().filter(x -> x.fragmentCoordinates().PositionUpper != NO_POSITION).collect(Collectors.toList());
        allUmiGroups.clear();
        allUmiGroups.addAll(filteredUmiGroups);

        Set<FragmentCoords> forwardCoordSet = Sets.newHashSet();
        Map<FragmentCoords, List<Integer>> groupIndicesByKey = Maps.newHashMap();
        UnionFind<Integer> mergedGroupIndices = new UnionFind<>();
        for(int i = 0; i < allUmiGroups.size(); i++)
        {
            mergedGroupIndices.add(i);

            DuplicateGroup umiGroup = allUmiGroups.get(i);
            Set<FragmentCoords> originalFragmentCoordinates = umiGroup.preCollapsedFragmentCoordinates();
            for(FragmentCoords coord : originalFragmentCoordinates)
            {
                FragmentCoords forwardCoord = coord.withFragmentOrientation(FORWARD);
                forwardCoordSet.add(forwardCoord);
                groupIndicesByKey.computeIfAbsent(forwardCoord, x -> Lists.newArrayList());
                groupIndicesByKey.get(forwardCoord).add(i);
            }
        }

        List<FragmentCoords> uniqueForwardCoords = Lists.newArrayList(forwardCoordSet);
        UnionFind<Integer> mergedForwardCoordsIndices = new UnionFind<>();
        for(int i = 0; i < uniqueForwardCoords.size(); i++)
            mergedForwardCoordsIndices.add(i);

        for(int i = 0; i < uniqueForwardCoords.size() - 1; i++)
        {
            FragmentCoords firstCoord = uniqueForwardCoords.get(i);
            for(int j = i + 1; j < uniqueForwardCoords.size(); j++)
            {
                FragmentCoords secondCoord = uniqueForwardCoords.get(j);
                if((firstCoord.PositionUpper == 0) != (secondCoord.PositionUpper == 0))
                    continue;

                if(firstCoord.PositionUpper == 0)
                {
                    if(!firstCoord.ChromsomeLower.equals(secondCoord.ChromsomeLower))
                        continue;

                    if(firstCoord.OrientLower != secondCoord.OrientLower)
                        continue;

                    if((firstCoord.SuppReadInfo != null) != (secondCoord.SuppReadInfo != null))
                        continue;

                    if(firstCoord.ReadIsLower != secondCoord.ReadIsLower)
                        continue;

                    if(firstCoord.UnmappedSourced != secondCoord.UnmappedSourced)
                        continue;

                    int lowerDiff = abs(firstCoord.PositionLower - secondCoord.PositionLower);
                    if(!(lowerDiff <= 10))
                        continue;

                    mergedForwardCoordsIndices.merge(i, j);
                    continue;
                }

                if(!firstCoord.ChromsomeLower.equals(secondCoord.ChromsomeLower))
                    continue;

                if(!firstCoord.ChromsomeUpper.equals(secondCoord.ChromsomeUpper))
                    continue;

                if(firstCoord.OrientLower != secondCoord.OrientLower)
                    continue;

                if(firstCoord.OrientUpper != secondCoord.OrientUpper)
                    continue;

                if((firstCoord.SuppReadInfo != null) != (secondCoord.SuppReadInfo != null))
                    continue;

                if(firstCoord.ReadIsLower != secondCoord.ReadIsLower)
                    continue;

                int lowerDiff = abs(firstCoord.PositionLower - secondCoord.PositionLower);
                int upperDiff = abs(firstCoord.PositionUpper - secondCoord.PositionUpper);
                if(lowerDiff != 0 && upperDiff != 0)
                    continue;

                if(!(lowerDiff <= 10 && upperDiff <= 10))
                    continue;

                mergedForwardCoordsIndices.merge(i, j);
            }
        }

        for(Set<Integer> coordIndices : mergedForwardCoordsIndices.getPartitions())
        {
            if(coordIndices.size() == 1)
                continue;

            List<FragmentCoords> coords = coordIndices.stream().map(uniqueForwardCoords::get).collect(Collectors.toList());
            FragmentCoords firstCoord = coords.get(0);
            for(int j = 1; j < coords.size(); j++)
            {
                FragmentCoords secondCoord = coords.get(j);
                groupIndicesByKey.get(firstCoord).addAll(groupIndicesByKey.get(secondCoord));
                groupIndicesByKey.remove(secondCoord);
            }
        }

        for(List<Integer> indices : groupIndicesByKey.values())
        {
            for(int i = 1; i < indices.size(); i++)
            {
                mergedGroupIndices.merge(indices.get(0), indices.get(i));
            }
        }

        for(Set<Integer> mergedIndices : mergedGroupIndices.getPartitions())
        {
            List<DuplicateGroup> groups = mergedIndices.stream().sorted().map(allUmiGroups::get).collect(Collectors.toList());
            List<DuplicateGroup> forwardGroups = groups.stream().filter(x -> x.fragmentCoordinates().forwardFragment()).collect(Collectors.toList());
            List<DuplicateGroup> reverseGroups = groups.stream().filter(x -> !x.fragmentCoordinates().forwardFragment()).collect(Collectors.toList());

            if(reverseGroups.isEmpty())
            {
                finalUmiGroups.addAll(forwardGroups);
                continue;
            }

            if(forwardGroups.isEmpty())
            {
                finalUmiGroups.addAll(reverseGroups);
                continue;
            }

            for(DuplicateGroup first : forwardGroups)
            {
                DuplicateGroup firstGroup = first;
                String firstUmi = firstGroup.umiId();
                FragmentCoords firstFragCoords = firstGroup.fragmentCoordinates();

                int secondIndex = 0;
                while(secondIndex < reverseGroups.size())
                {
                    DuplicateGroup secondGroup = reverseGroups.get(secondIndex);
                    String secondUmi = secondGroup.umiId();

                    boolean canCollapse = mUmiConfig.Duplex ?
                            hasDuplexUmiMatch(firstUmi, secondUmi, mUmiConfig.DuplexDelim, mUmiConfig.PermittedBaseDiff) : false;

                    if(canCollapse)
                    {
                        // merge the two opposing fragments / groups
                        reverseGroups.remove(secondIndex);
                        firstGroup.addReads(secondGroup.reads());
                        firstGroup.registerDualStrand();

                        // collapsing only occurs between a pair, not 1:M
                        break;
                    }
                    else
                    {
                        ++secondIndex;
                    }
                }
            }

            finalUmiGroups.addAll(forwardGroups);
            finalUmiGroups.addAll(reverseGroups);
        }

        allUmiGroups.clear();
        allUmiGroups.addAll(finalUmiGroups);
    }

    @VisibleForTesting
    public static void collapsePolyGDuplexUmis(final SequencingType sequencingType, final UmiConfig umiConfig,
            final List<DuplicateGroup> umiGroups, final List<ReadInfo> singleFragments)
    {
        if(!(sequencingType == SequencingType.ILLUMINA || sequencingType == SequencingType.BIOMODAL))
            return;

        if(!umiConfig.Duplex)
            return;

        // early exit: there are no fragments with unmapped mates
        boolean hasUnmappedMateFragments = Stream.concat(
                        umiGroups.stream().map(DuplicateGroup::fragmentCoordinates),
                        singleFragments.stream().map(ReadInfo::coordinates))
                .anyMatch(x -> !x.UnmappedSourced && x.PositionUpper == NO_POSITION);
        if(!hasUnmappedMateFragments)
            return;

        // add single fragments as umiGroups
        for(ReadInfo readInfo : singleFragments)
        {
            SAMRecord read = readInfo.read();
            FragmentCoords currentCoords = readInfo.coordinates();
            FragmentCoords preCollapsedCoords = readInfo.preCollapsedCoordinates();
            String umiId = umiConfig.extractUmiId(read.getReadName());
            DuplicateGroup duplicateGroup = new DuplicateGroup(umiId, Lists.newArrayList(read), preCollapsedCoords);
            duplicateGroup.updateFragmentCoordinates(currentCoords);
            umiGroups.add(duplicateGroup);
        }

        singleFragments.clear();

        // find candidate umi groups with unmapped mate
        UnionFind<Integer> umiGroupMerger = new UnionFind<>();
        SortedMap<Integer, SortedSet<Integer>> unmappedMateGroupsIndices = Maps.newTreeMap();
        for(int i = 0; i < umiGroups.size(); i++)
        {
            DuplicateGroup duplicateGroup = umiGroups.get(i);
            umiGroupMerger.add(i);
            FragmentCoords coords = duplicateGroup.fragmentCoordinates();
            if(coords.UnmappedSourced || coords.PositionUpper != NO_POSITION)
                continue;

            String umiId = duplicateGroup.umiId();
            if(polyGTailLength(umiId) < MIN_POLYG_UMI_TAIL_LENGTH)
                continue;

            for(FragmentCoords preCollapsedCoords : duplicateGroup.preCollapsedFragmentCoordinates())
            {
                int position = preCollapsedCoords.ReadIsLower ? preCollapsedCoords.PositionLower : preCollapsedCoords.PositionUpper;
                unmappedMateGroupsIndices.computeIfAbsent(position, key -> Sets.newTreeSet());
                unmappedMateGroupsIndices.get(position).add(i);
            }
        }

        // merge other umi groups into the above umi groups with unmapped mate
        for(int i = 0; i < umiGroups.size(); i++)
        {
            DuplicateGroup duplicateGroup = umiGroups.get(i);
            if(duplicateGroup.fragmentCoordinates().UnmappedSourced)
                continue;

            String umiId = duplicateGroup.umiId();

            for(FragmentCoords preCollapsedCoords : duplicateGroup.preCollapsedFragmentCoordinates())
            {
                int position = preCollapsedCoords.ReadIsLower ? preCollapsedCoords.PositionLower : preCollapsedCoords.PositionUpper;
                SortedSet<Integer> unmappedMateGroupIndices = unmappedMateGroupsIndices.get(position);
                if(unmappedMateGroupIndices == null)
                    continue;

                for(int j : unmappedMateGroupIndices)
                {
                    if(i == j)
                        continue;

                    DuplicateGroup unmappedMateGroup = umiGroups.get(j);
                    String unmappedMateUmiId = unmappedMateGroup.umiId();
                    String unmappedMateUmiIdPrefix = trimPolyGTail(unmappedMateUmiId);
                    String umiIdPrefix = trimPolyGTail(umiId);
                    int trimmedLength = min(unmappedMateUmiIdPrefix.length(), umiIdPrefix.length());
                    unmappedMateUmiIdPrefix = unmappedMateUmiIdPrefix.substring(0, trimmedLength);
                    umiIdPrefix = umiIdPrefix.substring(0, trimmedLength);
                    if(!umiIdPrefix.equals(unmappedMateUmiIdPrefix))
                        continue;

                    boolean merge = false;
                    for(FragmentCoords unmappedMateCoords : unmappedMateGroup.preCollapsedFragmentCoordinates())
                    {
                        int unmappedMateReadPosition = unmappedMateCoords.ReadIsLower
                                ? unmappedMateCoords.PositionLower
                                : unmappedMateCoords.PositionUpper;

                        if(unmappedMateReadPosition != position)
                            continue;

                        merge = true;
                        break;
                    }

                    if(merge)
                        umiGroupMerger.merge(i, j);
                }
            }
        }

        // reconstitute umiGroups and singleFragments based on the merged groups
        Collection<Set<Integer>> partitions = umiGroupMerger.getPartitions();
        List<DuplicateGroup> mergedUmiGroups = Lists.newArrayList();
        for(Set<Integer> partition : partitions)
        {
            List<DuplicateGroup> partitionGroups = partition.stream().map(umiGroups::get).toList();
            if(partitionGroups.size() == 1)
            {
                mergedUmiGroups.add(partitionGroups.get(0));
                continue;
            }

            int firstFullyMappedGroupIdx = 0;
            for(int i = 0; i < partitionGroups.size(); i++)
            {
                DuplicateGroup group = partitionGroups.get(i);
                if(group.fragmentCoordinates().PositionUpper != NO_POSITION)
                {
                    firstFullyMappedGroupIdx = i;
                    break;
                }
            }

            DuplicateGroup baseGroup = partitionGroups.get(firstFullyMappedGroupIdx);
            for(int i = 0; i < partitionGroups.size(); i++)
            {
                if(i == firstFullyMappedGroupIdx)
                    continue;

                baseGroup.addReads(partitionGroups.get(i).reads());
            }

            mergedUmiGroups.add(baseGroup);
        }

        umiGroups.clear();
        for(DuplicateGroup group : mergedUmiGroups)
        {
            if(group.readCount() == 1)
            {
                singleFragments.add(new ReadInfo(group.reads().get(0), group.fragmentCoordinates()));
                continue;
            }

            umiGroups.add(group);
        }
    }

    @VisibleForTesting
    public static boolean hasDuplexUmiMatch(final String first, final String second, final String duplexDelim, int permittedDiff)
    {
        String[] umiParts1 = splitUmi(first, duplexDelim);
        String[] umiParts2 = splitUmi(second, duplexDelim);

        if(umiParts1.length != 2 || umiParts2.length != 2)
            return false;

        return !exceedsUmiIdDiff(umiParts1[0], umiParts2[1], permittedDiff) && !exceedsUmiIdDiff(umiParts1[1], umiParts2[0], permittedDiff);
    }

    private static String[] splitUmi(final String duplexUmi, final String duplexDelim)
    {
        int delimIndex = duplexUmi.indexOf(duplexDelim);
        return new String[] { duplexUmi.substring(0, delimIndex), duplexUmi.substring(delimIndex + 1) };
    }
}

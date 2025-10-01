package com.hartwig.hmftools.redux.umi;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.collect.MergeUtils.clusterMerger;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ILLUMINA;
import static com.hartwig.hmftools.redux.ReduxConstants.MAX_IMBALANCED_UMI_BASE_DIFF;
import static com.hartwig.hmftools.redux.ReduxConstants.MAX_IMBALANCED_UMI_COUNT;
import static com.hartwig.hmftools.redux.ReduxConstants.MAX_UMI_BASE_DIFF_JITTER_COLLAPSE;
import static com.hartwig.hmftools.redux.ReduxConstants.MIN_POLYG_UMI_TAIL_LENGTH;
import static com.hartwig.hmftools.redux.common.DuplicateGroupCollapser.SINGLE_END_JITTER_COLLAPSE_DISTANCE;
import static com.hartwig.hmftools.redux.common.DuplicateGroupCollapser.collapseToNonOrientedKeyWithoutCoordinates;
import static com.hartwig.hmftools.redux.consensus.TemplateReads.selectTemplateRead;
import static com.hartwig.hmftools.redux.umi.UmiUtils.exceedsUmiIdDiff;
import static com.hartwig.hmftools.redux.umi.UmiUtils.polyGTailLength;
import static com.hartwig.hmftools.redux.umi.UmiUtils.trimPolyGTail;

import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
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
import com.hartwig.hmftools.common.collect.Union;
import com.hartwig.hmftools.common.collect.UnionFind;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.redux.common.DuplicateGroup;
import com.hartwig.hmftools.redux.common.FragmentCoords;
import com.hartwig.hmftools.redux.common.ReadInfo;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.samtools.SAMRecord;

public class UmiGroupBuilder
{
    private final SequencingType mSequencing;
    private final UmiConfig mUmiConfig;
    private final UmiStatistics mStats;

    public UmiGroupBuilder(final SequencingType sequencing, final UmiConfig config, final UmiStatistics stats)
    {
        mSequencing = sequencing;
        mUmiConfig = config;
        mStats = stats;
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

        Map<String,CoordinateGroup> coordinateGroupMap = formCoordGroups ? Maps.newHashMap() : null;

        List<DuplicateGroup> allUmiGroups = Lists.newArrayList();

        for(DuplicateGroup duplicateGroup : duplicateGroups)
        {
            List<DuplicateGroup> umiGroups = buildUmiGroups(duplicateGroup.fragmentCoordinates(), duplicateGroup.reads(), mUmiConfig);
            allUmiGroups.addAll(umiGroups);
        }

        collapsePolyGDuplexUmis(mSequencing, mUmiConfig, allUmiGroups, singleFragments);

        if(formCoordGroups)
        {
            // add in order of descending by fragment count for non-duplex collapsing
            Collections.sort(allUmiGroups, new UmiUtils.SizeComparator());
            for(DuplicateGroup umiGroup : allUmiGroups)
            {
                CoordinateGroup coordGroup = getOrCreateCoordGroup(coordinateGroupMap, umiGroup.fragmentCoordinates().keyNonOriented());
                coordGroup.addGroup(umiGroup);
            }

            allUmiGroups.clear();

            // add in single fragments
            for(ReadInfo readInfo : singleFragments)
            {
                CoordinateGroup coordGroup = getOrCreateCoordGroup(coordinateGroupMap, readInfo.coordinates().keyNonOriented());
                coordGroup.addSingleRead(readInfo);
            }

            if(mUmiConfig.BaseStats)
            {
                // test UMI similarity for all fragments and groups with the same coordinates
                for(CoordinateGroup coordGroup : coordinateGroupMap.values())
                {
                    captureUmiGroupStats(coordGroup.ForwardGroups);
                    captureUmiGroupStats(coordGroup.ReverseGroups);
                }
            }

            // collapse duplex and single UMIs with opposite orientations
            for(CoordinateGroup coordGroup : coordinateGroupMap.values())
            {
                collapseCoordinateGroup(coordGroup, allUmiGroups, singleFragments);
            }
        }

        jitterCollapseUmiGroups(mSequencing, mUmiConfig, allUmiGroups, singleFragments);

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

        if(captureStats)
        {
            Collection<CoordinateGroup> coordinateGroups = coordinateGroupMap != null ? coordinateGroupMap.values() : Collections.emptyList();
            captureStats(duplicateGroups, singleFragments, finalUmiGroups, formCoordGroups, coordinateGroups);
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
        List<DuplicateGroup> orderedGroups = groups.values()
                .stream()
                .sorted((new UmiUtils.SizeComparator()).thenComparing(DuplicateGroup::umiId))
                .collect(Collectors.toList());

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

    private class CoordinateGroup
    {
        public final String CoordKey;

        // store any mix of duplicate groups or single fragments
        public List<Object> ForwardGroups;
        public List<Object> ReverseGroups;

        public CoordinateGroup(final String coordKey)
        {
            CoordKey = coordKey;
            ForwardGroups = null;
            ReverseGroups = null;
        }

        public boolean hasOpposites()
        {
            return ForwardGroups != null && ReverseGroups != null;
        }

        private void addFragmentGroup(final Object object, boolean isForward)
        {
            if(isForward)
            {
                if(ForwardGroups == null)
                    ForwardGroups = Lists.newArrayList(object);
                else
                    ForwardGroups.add(object);
            }
            else
            {
                if(ReverseGroups == null)
                    ReverseGroups = Lists.newArrayList(object);
                else
                    ReverseGroups.add(object);
            }
        }

        public void addGroup(final DuplicateGroup group)
        {
            addFragmentGroup(group, group.fragmentCoordinates().forwardFragment());
        }

        public void addSingleRead(final ReadInfo readInfo)
        {
            addFragmentGroup(readInfo, readInfo.coordinates().forwardFragment());
        }

        public String toString()
        {
            return format("%s fwd(%d) rev(%d)",
                CoordKey, ForwardGroups != null ? ForwardGroups.size() : 0, ReverseGroups != null ? ReverseGroups.size() : 0);
        }
    }

    private CoordinateGroup getOrCreateCoordGroup(final Map<String,CoordinateGroup> coordinateGroups, final String coordKey)
    {
        CoordinateGroup coordinateGroup = coordinateGroups.get(coordKey);

        if(coordinateGroup != null)
            return coordinateGroup;

        CoordinateGroup newGroup = new CoordinateGroup(coordKey);
        coordinateGroups.put(coordKey, newGroup);
        return newGroup;
    }

    private void addUmiGroup(final List<DuplicateGroup> allUmiGroups, final List<Object> fragGroups)
    {
        if(fragGroups == null)
            return;

        for(Object fragGroup : fragGroups)
        {
            if(fragGroup instanceof DuplicateGroup)
            {
                allUmiGroups.add((DuplicateGroup) fragGroup);
            }
        }
    }

    private void collapseCoordinateGroup(
            final CoordinateGroup coordGroup, final List<DuplicateGroup> allUmiGroups, final List<ReadInfo> singleFragments)
    {
        // up until now fragments with the same coordinates but different orientation (ie F1R2 vs F2R1) have been kept separate,
        // but now merge these if they either don't use UMIs or have complementary duplex UMIs
        if(!coordGroup.hasOpposites() || !mUmiConfig.Duplex)
        {
            addUmiGroup(allUmiGroups, coordGroup.ForwardGroups);
            addUmiGroup(allUmiGroups, coordGroup.ReverseGroups);
            return;
        }

        for(Object first : coordGroup.ForwardGroups)
        {
            DuplicateGroup firstGroup = null;
            ReadInfo firstSingleRead = null;
            String firstUmi;
            FragmentCoords firstFragCoords;

            if(first instanceof DuplicateGroup)
            {
                firstGroup = (DuplicateGroup) first;
                firstUmi = firstGroup.umiId();
                firstFragCoords = firstGroup.fragmentCoordinates();
            }
            else
            {
                firstSingleRead = (ReadInfo) first;
                firstUmi = mUmiConfig.extractUmiId(firstSingleRead.id());
                firstFragCoords = firstSingleRead.coordinates();
            }

            int secondIndex = 0;
            while(secondIndex < coordGroup.ReverseGroups.size())
            {
                Object second = coordGroup.ReverseGroups.get(secondIndex);
                DuplicateGroup secondGroup = null;
                ReadInfo secondSingleRead = null;
                String secondUmi;

                if(second instanceof DuplicateGroup)
                {
                    secondGroup = (DuplicateGroup) second;
                    secondUmi = secondGroup.umiId();
                }
                else
                {
                    secondSingleRead = (ReadInfo) second;
                    secondUmi = mUmiConfig.extractUmiId(secondSingleRead.read().getReadName());
                }

                boolean canCollapse = mUmiConfig.Duplex ?
                        hasDuplexUmiMatch(firstUmi, secondUmi, mUmiConfig.DuplexDelim, mUmiConfig.PermittedBaseDiff) : false;

                if(canCollapse)
                {
                    // merge the two opposing fragments / groups
                    coordGroup.ReverseGroups.remove(secondIndex);

                    if(firstGroup == null)
                    {
                        // turn fragment into group
                        firstGroup = new DuplicateGroup(firstUmi, firstSingleRead.read(), firstFragCoords);

                        // and remove single fragment from the full set of these
                        singleFragments.remove(firstSingleRead);
                    }

                    if(secondGroup != null)
                    {
                        firstGroup.addReads(secondGroup.reads());
                        firstGroup.addPolyGUmiReads(secondGroup.polyGUmiReads());
                    }
                    else
                    {
                        firstGroup.addRead(secondSingleRead.read());

                        singleFragments.remove(secondSingleRead);
                    }

                    firstGroup.registerDualStrand();

                    // collapsing only occurs between a pair, not 1:M
                    break;
                }
                else
                {
                    ++secondIndex;
                }
            }

            if(firstGroup != null)
                allUmiGroups.add(firstGroup);
        }

        for(Object fragGroup : coordGroup.ReverseGroups)
        {
            if(fragGroup instanceof DuplicateGroup)
                allUmiGroups.add((DuplicateGroup)fragGroup);
        }
    }

    private static class UnmappedPair
    {
        public DuplicateGroup MappedGroup = null;
        public DuplicateGroup UnmappedGroup = null;
    }

    @VisibleForTesting
    public static void collapsePolyGDuplexUmis(final SequencingType sequencingType, final UmiConfig umiConfig,
            final List<DuplicateGroup> umiGroups, final List<ReadInfo> singleFragments)
    {
        if(umiGroups.isEmpty() && singleFragments.isEmpty())
            return;

        if(sequencingType != ILLUMINA)
            return;

        if(!umiConfig.Duplex)
            return;

        boolean hasPrimaryUnmappedMateFragments = Stream.concat(
                        umiGroups.stream().map(DuplicateGroup::fragmentCoordinates),
                        singleFragments.stream().map(ReadInfo::coordinates))
                .anyMatch(x -> !x.Unpaired && x.SuppReadInfo == null && x.PositionUpper == NO_POSITION);
        if(!hasPrimaryUnmappedMateFragments)
            return;

        // add non-supp single fragments as umiGroups
        List<ReadInfo> suppSingleFragments = Lists.newArrayList();
        for(ReadInfo readInfo : singleFragments)
        {
            FragmentCoords coords = readInfo.coordinates();
            if(coords.SuppReadInfo != null)
            {
                suppSingleFragments.add(readInfo);
                continue;
            }

            SAMRecord read = readInfo.read();
            String umiId = umiConfig.extractUmiId(read.getReadName());
            DuplicateGroup duplicateGroup = new DuplicateGroup(umiId, read, coords);
            umiGroups.add(duplicateGroup);
        }

        singleFragments.clear();
        singleFragments.addAll(suppSingleFragments);

        Map<String, UnmappedPair> pairedUnmappedGroups = Maps.newHashMap();
        List<DuplicateGroup> suppFragments = Lists.newArrayList();
        List<Union<DuplicateGroup, UnmappedPair>> fragments = Lists.newArrayList();
        UnionFind<Integer> umiGroupMerger = new UnionFind<>();
        for(DuplicateGroup group : umiGroups)
        {
            FragmentCoords coords = group.fragmentCoordinates();
            if(coords.SuppReadInfo != null)
            {
                suppFragments.add(group);
                continue;
            }

            if(coords.PositionUpper != NO_POSITION)
            {
                umiGroupMerger.add(fragments.size());
                fragments.add(Union.createLeft(group));
                continue;
            }

            String key = coords.lowerCoordinate() + ":" + group.umiId();
            pairedUnmappedGroups.computeIfAbsent(key, k -> new UnmappedPair());
            if(coords.UnmappedSourced)
            {
                pairedUnmappedGroups.get(key).UnmappedGroup = group;
                continue;
            }

            pairedUnmappedGroups.get(key).MappedGroup = group;
        }

        umiGroups.clear();
        umiGroups.addAll(suppFragments);

        HashMap<Pair<Integer, Boolean>, SortedSet<Integer>> groupedPolyGUnmappedPairIndices = Maps.newHashMap();
        for(UnmappedPair unmappedPair : pairedUnmappedGroups.values())
        {
            int i = fragments.size();
            fragments.add(Union.createRight(unmappedPair));
            umiGroupMerger.add(i);
            String umiId = unmappedPair.MappedGroup.umiId();
            if(polyGTailLength(umiId) < MIN_POLYG_UMI_TAIL_LENGTH)
                continue;

            FragmentCoords coords = unmappedPair.MappedGroup.fragmentCoordinates();
            int pos = coords.PositionLower;
            Orientation orient = coords.OrientLower;
            Pair<Integer, Boolean> key = Pair.of(pos, orient.isForward());
            groupedPolyGUnmappedPairIndices.computeIfAbsent(key, k -> Sets.newTreeSet());
            groupedPolyGUnmappedPairIndices.get(key).add(i);
        }

        // merge the above umi groups with unmapped mates into other umi groups
        for(int i = 0; i < fragments.size(); i++)
        {
            Union<DuplicateGroup, UnmappedPair> fragment = fragments.get(i);
            DuplicateGroup group = fragment.hasLeft() ? fragment.left() : fragment.right().MappedGroup;
            String umiId = group.umiId();
            FragmentCoords coords = group.fragmentCoordinates();
            int position = coords.PositionUpper == NO_POSITION || coords.ReadIsLower ? coords.PositionLower : coords.PositionUpper;
            Orientation orient = coords.PositionUpper == NO_POSITION || coords.ReadIsLower ? coords.OrientLower : coords.OrientUpper;
            Pair<Integer, Boolean> key = Pair.of(position, orient.isForward());
            SortedSet<Integer> polyGUnmappedPairIndices = groupedPolyGUnmappedPairIndices.get(key);
            if(polyGUnmappedPairIndices == null)
                continue;

            for(int j : polyGUnmappedPairIndices)
            {
                if(i == j)
                    continue;

                UnmappedPair unmappedPair = fragments.get(j).right();
                String unmappedUmiId = unmappedPair.MappedGroup.umiId();
                String unmappedUmiIdPrefix = trimPolyGTail(unmappedUmiId);
                String umiIdPrefix = trimPolyGTail(umiId);

                int trimmedLength = min(unmappedUmiIdPrefix.length(), umiIdPrefix.length());
                unmappedUmiIdPrefix = unmappedUmiIdPrefix.substring(0, trimmedLength);
                umiIdPrefix = umiIdPrefix.substring(0, trimmedLength);
                if(!umiIdPrefix.equals(unmappedUmiIdPrefix))
                    continue;

                umiGroupMerger.merge(i, j);
            }
        }

        // reconstitute umiGroups and singleFragments based on the merged groups
        Comparator<Union<DuplicateGroup, UnmappedPair>> baseFragmentComparator = Comparator
                .comparingInt((final Union<DuplicateGroup, UnmappedPair> x) -> x.hasLeft() ? 0 : 1)
                .thenComparing(Comparator.comparingInt((final Union<DuplicateGroup, UnmappedPair> x) ->
                        x.hasLeft()
                        ? x.left().readCount()
                        : x.right().MappedGroup.readCount()).reversed())
                .thenComparing((final Union<DuplicateGroup, UnmappedPair> x) -> selectTemplateRead(
                        x.hasLeft()
                        ? x.left().reads()
                        : x.right().MappedGroup.reads(),
                        x.hasLeft()
                        ? x.left().fragmentCoordinates()
                        : x.right().MappedGroup.fragmentCoordinates()).getReadName());

        Collection<Set<Integer>> partitions = umiGroupMerger.getPartitions();
        List<DuplicateGroup> mergedUmiGroups = Lists.newArrayList();
        for(Set<Integer> partition : partitions)
        {
            List<Union<DuplicateGroup, UnmappedPair>> partitionFragments = partition.stream().map(fragments::get).toList();
            if(partitionFragments.size() == 1)
            {
                Union<DuplicateGroup, UnmappedPair> fragment = partitionFragments.get(0);
                if(fragment.hasLeft())
                {
                    mergedUmiGroups.add(fragment.left());
                    continue;
                }

                UnmappedPair unmappedPair = fragment.right();

                if(unmappedPair.MappedGroup != null)
                    mergedUmiGroups.add(unmappedPair.MappedGroup);

                if(unmappedPair.UnmappedGroup != null)
                    mergedUmiGroups.add(unmappedPair.UnmappedGroup);

                continue;
            }

            Union<DuplicateGroup, UnmappedPair> baseFragment = Collections.min(partitionFragments, baseFragmentComparator);
            if(baseFragment.hasLeft())
            {
                mergedUmiGroups.add(baseFragment.left());
            }
            else
            {
                UnmappedPair unmappedPair = baseFragment.right();
                if(unmappedPair.MappedGroup != null)
                    mergedUmiGroups.add(unmappedPair.MappedGroup);

                if(unmappedPair.UnmappedGroup != null)
                    mergedUmiGroups.add(unmappedPair.UnmappedGroup);
            }

            for(int i = 0; i < partitionFragments.size(); i++)
            {
                Union<DuplicateGroup, UnmappedPair> fragment = partitionFragments.get(i);
                if(fragment == baseFragment)
                    continue;

                if(fragment.hasLeft())
                {
                    mergedUmiGroups.add(fragment.left());
                    continue;
                }

                if(baseFragment.hasLeft())
                {
                    DuplicateGroup baseGroup = baseFragment.left();
                    if(fragment.right().MappedGroup != null)
                        baseGroup.addPolyGUmiReads(fragment.right().MappedGroup.reads());

                    if(fragment.right().UnmappedGroup != null)
                        baseGroup.addPolyGUmiReads(fragment.right().UnmappedGroup.reads());
                }
                else
                {
                    UnmappedPair unmappedPair = baseFragment.right();
                    if(unmappedPair.MappedGroup != null && fragment.right().MappedGroup != null)
                        unmappedPair.MappedGroup.addPolyGUmiReads(fragment.right().MappedGroup.reads());

                    if(unmappedPair.UnmappedGroup != null && fragment.right().UnmappedGroup != null)
                        unmappedPair.UnmappedGroup.addPolyGUmiReads(fragment.right().UnmappedGroup.reads());
                }
            }
        }

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

    private static DuplicateGroup jitterMergeFn(final DuplicateGroup acc, final DuplicateGroup group)
    {
        acc.addNonConsensusReads(group.reads());
        acc.addNonConsensusReads(group.nonConsensusReads());
        acc.addPolyGUmiReads(group.polyGUmiReads());
        return acc;
    }

    private static final Comparator<DuplicateGroup> DUPLICATE_GROUP_COMPARATOR = Comparator.comparingInt(DuplicateGroup::readCount)
            .reversed()
            .thenComparing(x -> selectTemplateRead(x.reads(), x.fragmentCoordinates()).getReadName());

    private record JitterCollapseElementKey(long keyId, FragmentCoords coords, String umiId)
    {
        public boolean canMerge(final UmiConfig umiConfig, final JitterCollapseElementKey o)
        {
            int lowerDiff = abs(coords.PositionLower - o.coords.PositionLower);
            int upperDiff = abs(coords.PositionUpper - o.coords.PositionUpper);
            if(lowerDiff > 0 && upperDiff > 0)
            {
                return false;
            }

            if(lowerDiff > SINGLE_END_JITTER_COLLAPSE_DISTANCE || upperDiff > SINGLE_END_JITTER_COLLAPSE_DISTANCE)
            {
                return false;
            }

            if(coords.FragmentOrient == o.coords.FragmentOrient)
            {
                return !exceedsUmiIdDiff(umiId, o.umiId, MAX_UMI_BASE_DIFF_JITTER_COLLAPSE);
            }

            if(!umiConfig.Duplex)
            {
                return false;
            }

            return hasDuplexUmiMatch(umiId, o.umiId, umiConfig.DuplexDelim, MAX_UMI_BASE_DIFF_JITTER_COLLAPSE);
        }
    }

    private static void jitterCollapseUmiGroups(
            final SequencingType sequencingType, final UmiConfig umiConfig, final List<DuplicateGroup> umiGroups,
            final List<ReadInfo> singleFragments)
    {
        if(umiGroups.isEmpty() && singleFragments.isEmpty())
            return;

        if(sequencingType != ILLUMINA)
            return;

        boolean isUnpaired = Stream.concat(
                        umiGroups.stream().map(DuplicateGroup::fragmentCoordinates),
                        singleFragments.stream().map(ReadInfo::coordinates)).anyMatch(x -> x.Unpaired);

        if(isUnpaired)
            return;

        for(ReadInfo readInfo : singleFragments)
        {
            SAMRecord read = readInfo.read();
            FragmentCoords coords = readInfo.coordinates();
            String umiId = umiConfig.extractUmiId(read.getReadName());
            DuplicateGroup duplicateGroup = new DuplicateGroup(umiId, read, coords);
            umiGroups.add(duplicateGroup);
        }

        singleFragments.clear();

        List<DuplicateGroup> finalUmiGroups = Lists.newArrayList();
        Map<String, List<DuplicateGroup>> keyGroups = Maps.newHashMap();
        for(DuplicateGroup group : umiGroups)
        {
            FragmentCoords coords = group.fragmentCoordinates();
            if(coords.PositionUpper == NO_POSITION)
            {
                finalUmiGroups.add(group);
                continue;
            }

            String collapsedKey = collapseToNonOrientedKeyWithoutCoordinates(coords);
            keyGroups.computeIfAbsent(collapsedKey, key -> Lists.newArrayList());
            keyGroups.get(collapsedKey).add(group);
        }

        umiGroups.clear();

        long keyId = 0;
        for(List<DuplicateGroup> keyGroup : keyGroups.values())
        {
            Map<JitterCollapseElementKey, DuplicateGroup> elements = Maps.newHashMap();
            SortedMap<Integer, SortedMap<Integer, HashSet<JitterCollapseElementKey>>> lowerToUpper = Maps.newTreeMap();
            SortedMap<Integer, SortedMap<Integer, HashSet<JitterCollapseElementKey>>> upperToLower = Maps.newTreeMap();
            for(DuplicateGroup group : keyGroup)
            {
                FragmentCoords coords = group.fragmentCoordinates();
                String umiId = group.umiId();
                JitterCollapseElementKey key = new JitterCollapseElementKey(keyId, coords, umiId);
                keyId++;
                elements.put(key, group);

                lowerToUpper.computeIfAbsent(coords.PositionLower, k -> Maps.newTreeMap());
                lowerToUpper.get(coords.PositionLower).computeIfAbsent(coords.PositionUpper, k -> Sets.newHashSet());
                lowerToUpper.get(coords.PositionLower).get(coords.PositionUpper).add(key);

                upperToLower.computeIfAbsent(coords.PositionUpper, k -> Maps.newTreeMap());
                upperToLower.get(coords.PositionUpper).computeIfAbsent(coords.PositionLower, k -> Sets.newHashSet());
                upperToLower.get(coords.PositionUpper).get(coords.PositionLower).add(key);
            }

            Map<JitterCollapseElementKey, Set<JitterCollapseElementKey>> mergeGroups = Maps.newHashMap();
            for(JitterCollapseElementKey key1 : elements.keySet())
            {
                int lower = key1.coords.PositionLower;
                int upper = key1.coords.PositionUpper;

                SortedMap<Integer, HashSet<JitterCollapseElementKey>> inRange = lowerToUpper.get(lower)
                        .subMap(upper - SINGLE_END_JITTER_COLLAPSE_DISTANCE, upper + SINGLE_END_JITTER_COLLAPSE_DISTANCE + 1);
                for(JitterCollapseElementKey key2 : inRange.values().stream().flatMap(HashSet::stream).toList())
                {
                    if(key1.equals(key2))
                        continue;

                    mergeGroups.computeIfAbsent(key1, k -> Sets.newHashSet());
                    mergeGroups.get(key1).add(key2);
                }

                inRange = upperToLower.get(upper)
                        .subMap(lower - SINGLE_END_JITTER_COLLAPSE_DISTANCE, lower + SINGLE_END_JITTER_COLLAPSE_DISTANCE + 1);
                for(JitterCollapseElementKey key2 : inRange.values().stream().flatMap(HashSet::stream).toList())
                {
                    if(key1.equals(key2))
                        continue;

                    mergeGroups.computeIfAbsent(key1, k -> Sets.newHashSet());
                    mergeGroups.get(key1).add(key2);
                }
            }

            List<DuplicateGroup> mergedGroups = clusterMerger(
                    elements, (x, y) -> x.canMerge(umiConfig, y), DUPLICATE_GROUP_COMPARATOR, UmiGroupBuilder::jitterMergeFn, mergeGroups);
            finalUmiGroups.addAll(mergedGroups);
        }

        umiGroups.addAll(finalUmiGroups);
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

    private void captureStats(
            final List<DuplicateGroup> duplicateGroups, final List<ReadInfo> singleFragments, final List<DuplicateGroup> finalUmiGroups,
            boolean formCoordGroups, final Collection<CoordinateGroup> coordinateGroups)
    {
        int uniqueCoordCount = 0; // count of distinct coordinates after collapsing

        // count of final primary fragments (UMIs and singles)
        int uniqueFragmentCount = 0;
        int maxCoordUmiCount = 0;

        if(formCoordGroups)
        {
            for(CoordinateGroup coordGroup : coordinateGroups)
            {
                if(coordGroup.ForwardGroups != null)
                {
                    uniqueFragmentCount += coordGroup.ForwardGroups.size();
                    ++uniqueCoordCount;
                    maxCoordUmiCount = max(maxCoordUmiCount, coordGroup.ForwardGroups.size());
                }

                if(coordGroup.ReverseGroups != null)
                {
                    uniqueFragmentCount += coordGroup.ReverseGroups.size();

                    if(!coordGroup.ReverseGroups.isEmpty())
                        ++uniqueCoordCount;

                    maxCoordUmiCount = max(maxCoordUmiCount, coordGroup.ReverseGroups.size());
                }
            }
        }
        else
        {
            // should be 1 distinct coordinate
            int ungroupedFragmentCount = singleFragments.size();
            uniqueFragmentCount = ungroupedFragmentCount + finalUmiGroups.size();
            uniqueCoordCount = duplicateGroups.size() + ungroupedFragmentCount;
        }

        int maxUmiFragmentCount = 0;
        DuplicateGroup maxUmiGroup = null;

        for(DuplicateGroup umiGroup : finalUmiGroups)
        {
            ++mStats.UmiGroups;

            if(umiGroup.readCount() > maxUmiFragmentCount)
            {
                maxUmiGroup = umiGroup;
                maxUmiFragmentCount = umiGroup.readCount();
            }
        }

        if(uniqueCoordCount > 0 && uniqueFragmentCount > 0)
            mStats.recordFragmentPositions(uniqueCoordCount, uniqueFragmentCount, maxCoordUmiCount, maxUmiGroup);
    }

    private void captureUmiGroupStats(final List<Object> fragGroups)
    {
        if(fragGroups == null)
            return;

        List<DuplicateGroup> groups = Lists.newArrayList();

        for(Object fragGroup : fragGroups)
        {
            if(fragGroup instanceof DuplicateGroup)
            {
                groups.add((DuplicateGroup) fragGroup);
            }
            else
            {
                SAMRecord read = (SAMRecord)fragGroup;
                groups.add(new DuplicateGroup(mUmiConfig.extractUmiId(read.getReadName()), read, null));
            }
        }

        mStats.recordUmiBaseStats(mUmiConfig, groups);
    }
}

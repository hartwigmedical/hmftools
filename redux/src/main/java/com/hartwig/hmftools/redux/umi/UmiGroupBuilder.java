package com.hartwig.hmftools.redux.umi;

import static java.lang.Math.max;

import static com.hartwig.hmftools.redux.common.Constants.MAX_IMBALANCED_UMI_BASE_DIFF;
import static com.hartwig.hmftools.redux.common.Constants.MAX_IMBALANCED_UMI_COUNT;
import static com.hartwig.hmftools.redux.umi.UmiUtils.exceedsUmiIdDiff;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.redux.common.DuplicateGroup;
import com.hartwig.hmftools.redux.common.FragmentCoords;
import com.hartwig.hmftools.redux.common.FragmentStatus;
import com.hartwig.hmftools.redux.common.ReadInfo;

import htsjdk.samtools.SAMRecord;

public class UmiGroupBuilder
{
    private final UmiConfig mUmiConfig;
    private final UmiStatistics mStats;

    public UmiGroupBuilder(final UmiConfig config, final UmiStatistics stats)
    {
        mUmiConfig = config;
        mStats = stats;
    }

    public List<DuplicateGroup> processUmiGroups(
            final List<DuplicateGroup> duplicateGroups, final List<ReadInfo> singleFragments, boolean captureStats)
    {
        // organise groups by their UMIs, applying base-difference collapsing rules
        // UMI stats require evaluation of uncollapsed UMI groups with the same coordinates
        // at the same time organise UMI groups by the coordinates

        boolean formCoordGroups = mUmiConfig.BaseStats || (duplicateGroups.size() + singleFragments.size() > 1);

        List<CoordinateGroup> coordinateGroups = formCoordGroups ? Lists.newArrayList() : null;

        List<DuplicateGroup> allUmiGroups = Lists.newArrayList();

        for(DuplicateGroup duplicateGroup : duplicateGroups)
        {
            List<DuplicateGroup> umiGroups = buildUmiGroups(duplicateGroup.fragmentCoordinates(), duplicateGroup.reads(), mUmiConfig);

            if(formCoordGroups)
            {
                CoordinateGroup coordGroup = getOrCreateCoordGroup(coordinateGroups, duplicateGroup.fragmentCoordinates().Key);

                // add in order of descending by fragment count for non-duplex collapsing
                Collections.sort(umiGroups, new UmiUtils.SizeComparator());
                umiGroups.forEach(x -> coordGroup.addGroup(x));
            }
            else
            {
                allUmiGroups.addAll(umiGroups);
            }
        }

        if(formCoordGroups)
        {
            // add in single fragments
            for(ReadInfo readInfo : singleFragments)
            {
                CoordinateGroup coordGroup = getOrCreateCoordGroup(coordinateGroups, readInfo.coordinates().Key);
                coordGroup.addSingleRead(readInfo);
            }

            if(mUmiConfig.BaseStats)
            {
                // test UMI similarity for all fragments and groups with the same coordinates
                for(CoordinateGroup coordGroup : coordinateGroups)
                {
                    captureUmiGroupStats(coordGroup.ForwardGroups);
                    captureUmiGroupStats(coordGroup.ReverseGroups);
                }
            }

            // collapse duplex and single UMIs with opposite orientations
            for(CoordinateGroup coordGroup : coordinateGroups)
            {
                collapseCoordinateGroup(allUmiGroups, coordGroup);
            }
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

        // remove singles which were collapsed into a UMI group
        int index = 0;
        while(index < singleFragments.size())
        {
            SAMRecord read = singleFragments.get(index).read();

            if(finalUmiGroups.stream().anyMatch(x -> x.reads().contains(read)))
                singleFragments.remove(index);
            else
                ++index;
        }

        if(captureStats)
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
    }

    private CoordinateGroup getOrCreateCoordGroup(final List<CoordinateGroup> coordinateGroups, final String coordKey)
    {
        for(CoordinateGroup coordinateGroup : coordinateGroups)
        {
            if(coordinateGroup.CoordKey.equals(coordKey))
                return coordinateGroup;
        }

        CoordinateGroup newGroup = new CoordinateGroup(coordKey);
        coordinateGroups.add(newGroup);
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

    private void collapseCoordinateGroup(final List<DuplicateGroup> allUmiGroups, final CoordinateGroup coordGroup)
    {
        // up until now fragments with the same coordinates but different ordering (ie F1R2 vs F2R1) have been kept separate.
        // but now merge these if they either don't use UMIs or have complementary duplex UMIs
        if(!coordGroup.hasOpposites())
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

                    if(firstGroup == null) // turn fragment into group
                    {
                        firstGroup = new DuplicateGroup(firstUmi, firstSingleRead.read(), firstFragCoords);
                    }

                    if(secondGroup != null)
                    {
                        firstGroup.addReads(secondGroup.reads());
                    }
                    else
                    {
                        firstGroup.addRead(secondSingleRead.read());
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

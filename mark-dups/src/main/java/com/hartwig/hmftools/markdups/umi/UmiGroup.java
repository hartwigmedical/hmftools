package com.hartwig.hmftools.markdups.umi;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UMI_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.orientation;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.Constants.MAX_UMI_BASE_DIFF;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.getUnclippedPosition;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.readToString;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.markdups.common.Fragment;
import com.hartwig.hmftools.markdups.common.FragmentCoordinates;

import htsjdk.samtools.SAMRecord;

public class UmiGroup
{
    private final String mId;
    private final List<Fragment> mFragments;
    private List<String> mReadIds;
    private int mFragmentCount;

    // reads from each fragment are organised into their like-types from which consensus reads can be formed
    private final List<SAMRecord>[] mReadGroups;
    private final boolean[] mReadGroupComplete;
    private final ReadTypeId[] mReadTypeIndex; // list of distinct read types
    private final FragmentCoordinates mFragmentCoordinates;

    private static final int MAX_READ_TYPES = 4;

    public UmiGroup(final String id, final Fragment fragment)
    {
        mId = id;
        mFragments = Lists.newArrayList(fragment);
        mReadIds = null;
        mReadGroups = new List[MAX_READ_TYPES];
        mReadGroupComplete = new boolean[MAX_READ_TYPES];
        mReadTypeIndex = new ReadTypeId[MAX_READ_TYPES];
        mFragmentCount = 0;
        mFragmentCoordinates = fragment.coordinates();
    }

    public List<Fragment> fragments() { return mFragments; }
    public int fragmentCount() { return mFragmentCount > 0 ? mFragmentCount : mFragments.size(); }
    public FragmentCoordinates fragmentCoordinates() { return mFragmentCoordinates; }

    public String id() { return mId; }

    public void categoriseReads()
    {
        if(mFragmentCount > 0)
            return;

        mFragmentCount = mFragments.size();
        mReadIds = Lists.newArrayListWithExpectedSize(mFragmentCount);

        // initialise lists
        Fragment firstFragment = mFragments.get(0);

        mReadGroups[ReadLegType.PRIMARY.ordinal()] = Lists.newArrayListWithExpectedSize(mFragmentCount);

        for(int i = 0; i < firstFragment.reads().size(); ++i)
        {
            SAMRecord read = firstFragment.reads().get(i);

            if(i == 0)
            {
                if(read.getReadPairedFlag())
                    mReadGroups[ReadLegType.MATE.ordinal()] = Lists.newArrayListWithExpectedSize(mFragmentCount);

                if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
                    mReadGroups[ReadLegType.PRIMARY_SUPPLEMENTARY.ordinal()] = Lists.newArrayListWithExpectedSize(mFragmentCount);
            }
            else if(!read.getSupplementaryAlignmentFlag() && read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
            {
                mReadGroups[ReadLegType.MATE_SUPPLEMENTARY.ordinal()] = Lists.newArrayListWithExpectedSize(mFragmentCount);
            }
        }

        for(Fragment fragment : mFragments)
        {
            mReadIds.add(fragment.id());
            fragment.setUmi(mId);
            fragment.reads().forEach(x -> addRead(x));
            fragment.reads().forEach(x -> x.setAttribute(UMI_ATTRIBUTE, id()));
        }

        mFragments.clear();
    }

    public void addRead(final SAMRecord read)
    {
        int readTypeIndex = getReadTypeIndex(read);

        if(mReadGroups[readTypeIndex] == null)
        {
            mReadGroups[readTypeIndex] = Lists.newArrayListWithExpectedSize(mFragmentCount);
        }

        mReadGroups[readTypeIndex].add(read);
    }

    private class ReadTypeId
    {
        public final String Chromosome;
        public final int Position; // unclipped pos if not supplementary
        public final byte Orientation;
        public final boolean FirstInPair;
        public final boolean Supplementary;
        public final boolean Unmapped;

        public ReadTypeId(
                final String chromosome, final int position, final byte orientation, final boolean supplementary,
                final boolean firstInPair, final boolean unmapped)
        {
            Chromosome = chromosome;
            Position = position;
            Orientation = orientation;
            Supplementary = supplementary;
            FirstInPair = firstInPair;
            Unmapped = unmapped;
        }

        public boolean matches(final SAMRecord read)
        {
            if(Unmapped || read.getReadUnmappedFlag())
                return Unmapped == read.getReadUnmappedFlag();

            if(!read.getReferenceName().equals(Chromosome))
                return false;

            if(orientation(read) != Orientation)
                return false;

            if(read.getSupplementaryAlignmentFlag() != Supplementary)
                return false;

            int unclippedPosition = getUnclippedPosition(read);

            if(Supplementary)
            {
                return abs(unclippedPosition - Position) < read.getBaseQualities().length;
            }
            else
            {
                return unclippedPosition == Position;
            }
        }

        public String toString()
        {
            return format("%s:%d%s%s %s",
                Chromosome, Position, Orientation == POS_ORIENT ? " " : "_R ",
                    FirstInPair ? "R1" : "R2", Supplementary ? "supp" : "");
        }
    }

    private int getReadTypeIndex(final SAMRecord read)
    {
        int nextIndex = 0;
        for(int i = 0; i < mReadTypeIndex.length; ++i)
        {
            if(mReadTypeIndex[i] == null)
            {
                nextIndex = i;
                break;
            }

            if(mReadTypeIndex[i].matches(read))
            {
                if(mReadTypeIndex[i].FirstInPair != read.getFirstOfPairFlag())
                {
                    MD_LOGGER.warn("umiGroup({}) read({}) readTypeIndex({}) mismatch",
                            mId, readToString(read), mReadTypeIndex[i]);
                }

                return i;
            }
        }

        if(mReadTypeIndex[MAX_READ_TYPES-1] != null)
        {
            // shouldn't happen but revert to matching purely on flag attributes
            for(int i = 0; i < mReadTypeIndex.length; ++i)
            {
                if(mReadTypeIndex[i].FirstInPair == read.getFirstOfPairFlag()
                && mReadTypeIndex[i].Supplementary == read.getSupplementaryAlignmentFlag())
                {
                    return i;
                }
            }

            return 0;
        }

        mReadTypeIndex[nextIndex] = new ReadTypeId(
                read.getReferenceName(),
                read.getReadUnmappedFlag() ? 0 : getUnclippedPosition(read),
                read.getReadUnmappedFlag() ? 0 : orientation(read),
                read.getSupplementaryAlignmentFlag(), read.getFirstOfPairFlag(), read.getReadUnmappedFlag());

        return nextIndex;
    }

    public boolean allReadsReceived()
    {
        for(int i = 0; i < mReadGroups.length; ++i)
        {
            if(mReadGroupComplete[i])
                continue;

            List<SAMRecord> readGroup = mReadGroups[i];

            if(readGroup != null && readGroup.size() < mFragmentCount)
                return false;
        }

        return true;
    }

    public boolean hasCompleteReadGroup()
    {
        return Arrays.stream(mReadGroups).anyMatch(x -> x != null && !x.isEmpty() && x.size() == mFragmentCount);
    }

    public List<SAMRecord> popCompletedReads(final ConsensusReads consensusReads, boolean processIncompletes)
    {
        List<SAMRecord> reads = Lists.newArrayList();

        for(int i = 0; i < mReadGroups.length; ++i)
        {
            if(mReadGroupComplete[i])
                continue;

            List<SAMRecord> readGroup = mReadGroups[i];

            if(readGroup == null || readGroup.isEmpty())
                continue;

            if(readGroup.size() < mFragmentCount && !processIncompletes)
                continue;

            reads.addAll(readGroup);

            ConsensusReadInfo consensusReadInfo;
            if(i == ReadLegType.PRIMARY_SUPPLEMENTARY.ordinal() || i == ReadLegType.MATE_SUPPLEMENTARY.ordinal())
            {
                consensusReadInfo = consensusReads.createConsensusRead(findConsistentSupplementaries(readGroup), mId);
            }
            else
            {
                consensusReadInfo = consensusReads.createConsensusRead(readGroup, mId);
            }

            reads.add(consensusReadInfo.ConsensusRead);

            readGroup.clear();
            mReadGroupComplete[i] = true;
        }

        return reads;
    }

    private List<SAMRecord> findConsistentSupplementaries(final List<SAMRecord> readGroup)
    {
        Map<Integer,Integer> posCounts = Maps.newHashMap();
        for(SAMRecord read : readGroup)
        {
            int unclippedPos = getUnclippedPosition(read);
            Integer count = posCounts.get(unclippedPos);
            posCounts.put(unclippedPos, count != null ? count + 1 : 1);
        }

        if(posCounts.size() == 1)
            return readGroup;

        int maxPos = 0;
        int maxCount = 0;
        for(Map.Entry<Integer,Integer> entry : posCounts.entrySet())
        {
            if(entry.getValue() > maxCount)
            {
                maxPos = entry.getKey();
                maxCount = entry.getValue();
            }
        }

        int finalMaxPos = maxPos;
        return readGroup.stream().filter(x -> getUnclippedPosition(x) == finalMaxPos).collect(Collectors.toList());
    }

    public List<String> getReadIds()
    {
        if(!mFragments.isEmpty())
            return mFragments.stream().map(x -> x.id()).collect(Collectors.toList());
        else
            return mReadIds;
    }

    private enum ReadLegType
    {
        PRIMARY,
        MATE,
        PRIMARY_SUPPLEMENTARY,
        MATE_SUPPLEMENTARY;
    }

    public int cachedReadCount()
    {
        return Arrays.stream(mReadGroups).filter(x -> x != null).mapToInt(x -> x.size()).sum();
    }

    public String toString()
    {
        if(mFragmentCount == 0)
            return format("id(%s) fragments(%d)", mId, mFragments.size());

        StringJoiner sj = new StringJoiner(", ");
        for(ReadLegType legType : ReadLegType.values())
        {
            List<SAMRecord> readGroup = mReadGroups[legType.ordinal()];

            if(readGroup == null)
                continue;

            sj.add(format("%s=%d", legType, mReadGroupComplete[legType.ordinal()] ? mFragmentCount : readGroup.size()));
        }

        return format("id(%s) fragments(%d) readCounts(%s)", mId, mFragmentCount, sj);
    }

    public static List<UmiGroup> buildUmiGroups(final List<Fragment> fragments, final UmiConfig config)
    {
        Map<String,UmiGroup> groups = Maps.newHashMap();
        boolean checkDefinedUmis = config.hasDefinedUmis();
        boolean useDefinedUmis = checkDefinedUmis;

        for(Fragment fragment : fragments)
        {
            String umiId = config.extractUmiId(fragment.id());

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

            UmiGroup group = groups.get(umiId);

            if(group == null)
            {
                groups.put(umiId, new UmiGroup(umiId, fragment));
            }
            else
            {
                group.fragments().add(fragment);
            }
        }

        if(useDefinedUmis)
        {
            return groups.values().stream().collect(Collectors.toList());
        }

        List<UmiGroup> orderedGroups = groups.values().stream().sorted(new SizeComparator()).collect(Collectors.toList());
        List<UmiGroup> finalGroups = Lists.newArrayList();

        int i = 0;
        while(i < orderedGroups.size() - 1)
        {
            UmiGroup first = orderedGroups.get(i);

            List<UmiGroup> cluster = Lists.newArrayList(first);

            int j = i + 1;
            while(j < orderedGroups.size())
            {
                UmiGroup second = orderedGroups.get(j);

                boolean merged = false;

                for(UmiGroup existing : cluster)
                {
                    if(existing.fragmentCount() >= second.fragmentCount() && !exceedsUmiIdDiff(existing.id(), second.id()))
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
                first.fragments().addAll(cluster.get(j).fragments());
            }

            finalGroups.add(first);
            ++i;
        }

        return orderedGroups;
    }

    public static boolean exceedsUmiIdDiff(final String first, final String second)
    {
        if(first.length() != second.length())
            return true;

        short diffs = 0;
        for(short i = 0; i < first.length(); ++i)
        {
            if(first.charAt(i) != second.charAt(i))
            {
                ++diffs;

                if(diffs > MAX_UMI_BASE_DIFF)
                    return true;
            }
        }

        return false;
    }

    private static class SizeComparator implements Comparator<UmiGroup>
    {
        public int compare(final UmiGroup first, final UmiGroup second)
        {
            if(first.fragments().size() < second.fragments().size())
                return 1;
            else if(first.fragments().size() > second.fragments().size())
                return -1;
            else
                return 0;
        }
    }
}

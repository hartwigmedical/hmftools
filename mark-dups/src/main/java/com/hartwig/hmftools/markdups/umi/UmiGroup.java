package com.hartwig.hmftools.markdups.umi;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UMI_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.orientation;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
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
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.markdups.common.Fragment;

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
    private final ReadTypeId[] mPrimaryReadTypeIndex; // details for primary and mate reads
    private final String mCoordinatesKey;

    private static final int MAX_READ_TYPES = ReadType.values().length;
    private static final int PRIMARY_READ_TYPES = ReadType.MATE.ordinal() + 1;

    public UmiGroup(final String id, final Fragment fragment)
    {
        mId = id;
        mFragments = Lists.newArrayList(fragment);
        mReadIds = null;
        mReadGroups = new List[MAX_READ_TYPES];
        mReadGroupComplete = new boolean[MAX_READ_TYPES];
        mPrimaryReadTypeIndex = new ReadTypeId[PRIMARY_READ_TYPES];
        mFragmentCount = 0;
        mCoordinatesKey = fragment.coordinates().Key;
    }

    public List<Fragment> fragments() { return mFragments; }
    public int fragmentCount() { return mFragmentCount > 0 ? mFragmentCount : mFragments.size(); }
    public String coordinatesKey() { return mCoordinatesKey; }

    public String id() { return mId; }

    public void categoriseReads()
    {
        if(mFragmentCount > 0)
            return;

        mFragmentCount = mFragments.size();
        mReadIds = Lists.newArrayListWithExpectedSize(mFragmentCount);

        // establish expected lists by read type
        Fragment firstFragment = mFragments.get(0);

        mReadGroups[ReadType.PRIMARY.ordinal()] = Lists.newArrayListWithExpectedSize(mFragmentCount);

        for(int i = 0; i < firstFragment.reads().size(); ++i)
        {
            SAMRecord read = firstFragment.reads().get(i);

            SupplementaryReadData suppData = SupplementaryReadData.from(read);
            boolean hasValidSupp = suppData != null && HumanChromosome.contains(suppData.Chromosome);

            if(i == 0)
            {
                if(read.getReadPairedFlag() && HumanChromosome.contains(read.getMateReferenceName()))
                    mReadGroups[ReadType.MATE.ordinal()] = Lists.newArrayListWithExpectedSize(mFragmentCount);

                if(hasValidSupp)
                    mReadGroups[ReadType.PRIMARY_SUPPLEMENTARY.ordinal()] = Lists.newArrayListWithExpectedSize(mFragmentCount);
            }
            else if(!read.getSupplementaryAlignmentFlag() && hasValidSupp)
            {
                mReadGroups[ReadType.MATE_SUPPLEMENTARY.ordinal()] = Lists.newArrayListWithExpectedSize(mFragmentCount);
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

    private enum ReadType
    {
        PRIMARY,
        MATE,
        PRIMARY_SUPPLEMENTARY,
        MATE_SUPPLEMENTARY;
    }

    private class ReadTypeId
    {
        public final String Chromosome;
        public final int UnclippedPosition;
        public final int[] PositionRange;
        public final byte Orientation;
        public final boolean FirstInPair;
        public final boolean HasSupplementary;
        public final boolean Unmapped;

        public ReadTypeId(
                final String chromosome, final int unclippedPosition, final int position, final byte orientation,
                final boolean hasSupplementary, final boolean firstInPair, final boolean unmapped)
        {
            Chromosome = chromosome;
            UnclippedPosition = unclippedPosition;
            PositionRange = new int[] { position, position };
            Orientation = orientation;
            HasSupplementary = hasSupplementary;
            FirstInPair = firstInPair;
            Unmapped = unmapped;
        }

        public void updatePosition(final int position)
        {
            PositionRange[SE_START] = min(PositionRange[SE_START], position);
            PositionRange[SE_END] = max(PositionRange[SE_END], position);
        }

        public boolean primaryMatches(final SAMRecord read)
        {
            if(Unmapped || read.getReadUnmappedFlag())
                return Unmapped == read.getReadUnmappedFlag();

            if(!read.getReferenceName().equals(Chromosome))
                return false;

            if(orientation(read) != Orientation)
                return false;

            return getUnclippedPosition(read) == UnclippedPosition;
        }

        public boolean supplementaryMatches(final SAMRecord read, final SupplementaryReadData suppData)
        {
            if(!suppData.Chromosome.equals(Chromosome))
                return false;

            if(orientation(read) != Orientation)
                return false;

            return positionWithin(suppData.Position, PositionRange[SE_START], PositionRange[SE_END]);
        }

        public String toString()
        {
            return format("%s:%d%s%s %d-%d %s",
                Chromosome, UnclippedPosition, Orientation == POS_ORIENT ? " " : "_R ",
                    FirstInPair ? "R1" : "R2", PositionRange[SE_START], PositionRange[SE_END],
                    HasSupplementary ? "has-supp" : "");
        }
    }

    private int getReadTypeIndex(final SAMRecord read)
    {
        if(!read.getSupplementaryAlignmentFlag())
        {
            int index = 0;
            for(; index < mPrimaryReadTypeIndex.length; ++index)
            {
                if(mPrimaryReadTypeIndex[index] == null)
                    break;

                if(mPrimaryReadTypeIndex[index].primaryMatches(read))
                {
                    mPrimaryReadTypeIndex[index].updatePosition(read.getAlignmentStart());

                    /* rare and unimportant
                    if(mReadTypeIndex[i].FirstInPair != read.getFirstOfPairFlag())
                    {
                        MD_LOGGER.trace("umiGroup({}) read({} {}) readTypeIndex({}) mismatch",
                                mId, readToString(read), read.getFirstOfPairFlag() ? "R1" : "R2", mReadTypeIndex[i]);
                    }
                    */

                    return index;
                }
            }

            boolean hasSupplementary = read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE);

            mPrimaryReadTypeIndex[index] = new ReadTypeId(
                    read.getReferenceName(),
                    read.getReadUnmappedFlag() ? 0 : getUnclippedPosition(read), read.getAlignmentStart(),
                    read.getReadUnmappedFlag() ? 0 : orientation(read),
                    hasSupplementary, read.getFirstOfPairFlag(), read.getReadUnmappedFlag());

            // register the expected supplementary read group
            if(hasSupplementary)
            {
                int suppReadTypeIndex = index == ReadType.PRIMARY.ordinal() ?
                        ReadType.PRIMARY_SUPPLEMENTARY.ordinal() : ReadType.MATE_SUPPLEMENTARY.ordinal();

                mReadGroups[suppReadTypeIndex] = Lists.newArrayListWithExpectedSize(mFragmentCount);
            }

            return index;
        }

        boolean checkSuppData = Arrays.stream(mPrimaryReadTypeIndex).filter(x -> x != null && x.HasSupplementary).count() == 2;
        SupplementaryReadData suppData = checkSuppData ? SupplementaryReadData.from(read) : null;

        int index = 0;
        for(; index < mPrimaryReadTypeIndex.length; ++index)
        {
            if(mPrimaryReadTypeIndex[index] == null)
                continue;

            if(checkSuppData && mPrimaryReadTypeIndex[index].supplementaryMatches(read, suppData))
                break;
            else if(!checkSuppData && mPrimaryReadTypeIndex[index].HasSupplementary)
                break;
        }

        return index == ReadType.PRIMARY.ordinal() ? ReadType.PRIMARY_SUPPLEMENTARY.ordinal() : ReadType.MATE_SUPPLEMENTARY.ordinal();
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
        return Arrays.stream(mReadGroups).anyMatch(x -> x != null && !x.isEmpty() && x.size() >= mFragmentCount);
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

            try
            {
                ConsensusReadInfo consensusReadInfo;

                if(i == ReadType.PRIMARY_SUPPLEMENTARY.ordinal() || i == ReadType.MATE_SUPPLEMENTARY.ordinal())
                {
                    // supplementaries can go to difference places and some reads have more than one, so go with the most frequent
                    consensusReadInfo = consensusReads.createConsensusRead(findConsistentSupplementaries(readGroup), mId);
                }
                else
                {
                    consensusReadInfo = consensusReads.createConsensusRead(readGroup, mId);
                }

                reads.add(consensusReadInfo.ConsensusRead);
            }
            catch(Exception e)
            {
                MD_LOGGER.error("error forming consensus: umi({}) coords({})", toString());

                for(SAMRecord read : readGroup)
                {
                    MD_LOGGER.error("read: {}", readToString(read));
                }

                e.printStackTrace();
            }

            readGroup.clear();
            mReadGroupComplete[i] = true;
        }

        return reads;
    }

    private List<SAMRecord> findConsistentSupplementaries(final List<SAMRecord> readGroup)
    {
        Map<String,Integer> posCounts = Maps.newHashMap();

        String[] readChrPositions = new String[readGroup.size()];

        for(int i = 0; i < readGroup.size(); ++i)
        {
            SAMRecord read = readGroup.get(i);
            String chrPosStr = format("%s_%d", read.getReferenceName(), getUnclippedPosition(read));
            readChrPositions[i] = chrPosStr;
            Integer count = posCounts.get(chrPosStr);
            posCounts.put(chrPosStr, count != null ? count + 1 : 1);
        }

        if(posCounts.size() == 1)
            return readGroup;

        if(posCounts.size() == 2)
            return Lists.newArrayList(readGroup.get(0));

        String maxChrPos = "";
        int maxCount = 0;
        for(Map.Entry<String,Integer> entry : posCounts.entrySet())
        {
            if(entry.getValue() > maxCount)
            {
                maxChrPos = entry.getKey();
                maxCount = entry.getValue();
            }
        }

        List<SAMRecord> selectedReads = Lists.newArrayListWithExpectedSize(maxCount);
        for(int i = 0; i < readGroup.size(); ++i)
        {
            if(readChrPositions[i].equals(maxChrPos))
                selectedReads.add(readGroup.get(i));
        }

        return selectedReads;
    }

    public List<String> getReadIds()
    {
        if(!mFragments.isEmpty())
            return mFragments.stream().map(x -> x.id()).collect(Collectors.toList());
        else
            return mReadIds;
    }

    public int cachedReadCount()
    {
        return Arrays.stream(mReadGroups).filter(x -> x != null).mapToInt(x -> x.size()).sum();
    }

    public String toString()
    {
        if(mFragmentCount == 0)
            return format("id(%s) fragments(%d) coords(%s)", mId, mFragments.size(), mCoordinatesKey);

        StringJoiner sj = new StringJoiner(", ");
        for(ReadType readType : ReadType.values())
        {
            int readTypeIndex = readType.ordinal();
            List<SAMRecord> readGroup = mReadGroups[readTypeIndex];

            if(readGroup == null)
                continue;

            sj.add(format("%s=%d %s", readType, readGroup.size(), mReadGroupComplete[readTypeIndex] ? "complete" : "pending"));
        }

        return format("id(%s) fragments(%d) coords(%s) readCounts(%s)", mId, mFragmentCount, mCoordinatesKey, sj);
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

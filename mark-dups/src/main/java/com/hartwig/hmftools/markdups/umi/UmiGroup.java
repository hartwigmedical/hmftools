package com.hartwig.hmftools.markdups.umi;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UMI_ATTRIBUTE;
import static com.hartwig.hmftools.markdups.common.Constants.MAX_UMI_BASE_DIFF;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.getUnclippedPosition;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.markdups.common.Fragment;

import htsjdk.samtools.SAMRecord;

public class UmiGroup
{
    private final String mId;
    private final List<Fragment> mFragments;
    private int mFragmentCount;

    private final List<SAMRecord>[] mReadGroups;
    private final boolean[] mReadGroupComplete;

    public UmiGroup(final String id, final Fragment fragment)
    {
        mId = id;
        mFragments = Lists.newArrayList(fragment);
        mReadGroups = new List[ReadLegType.values().length];
        mReadGroupComplete = new boolean[ReadLegType.values().length];
        mFragmentCount = 0;
    }

    public List<Fragment> fragments() { return mFragments; }
    public int fragmentCount() { return mFragmentCount > 0 ? mFragmentCount : mFragments.size(); }

    public String id() { return mId; }

    public void categoriseReads()
    {
        mFragmentCount = mFragments.size();

        // initialise lists
        Fragment firstFragment = mFragments.get(0);

        mReadGroups[ReadLegType.PRIMARY.ordinal()] = Lists.newArrayList();

        for(int i = 0; i < firstFragment.reads().size(); ++i)
        {
            SAMRecord read = firstFragment.reads().get(i);

            if(i == 0)
            {
                if(read.getReadPairedFlag())
                    mReadGroups[ReadLegType.MATE.ordinal()] = Lists.newArrayList();

                if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
                    mReadGroups[ReadLegType.PRIMARY_SUPPLEMENTARY.ordinal()] = Lists.newArrayList();
            }
            else if(!read.getSupplementaryAlignmentFlag() && read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
            {
                mReadGroups[ReadLegType.MATE_SUPPLEMENTARY.ordinal()] = Lists.newArrayList();
            }
        }

        for(Fragment fragment : mFragments)
        {
            fragment.setUmi(mId);
            fragment.reads().forEach(x -> addRead(x));
            fragment.reads().forEach(x -> x.setAttribute(UMI_ATTRIBUTE, id()));
        }
    }

    public void addRead(final SAMRecord read)
    {
        ReadLegType legType = getLegType(read);
        int legIndex = legType.ordinal();

        if(mReadGroups[legIndex] == null)
        {
            mReadGroups[legIndex] = Lists.newArrayList(read);
        }
        else
        {
            mReadGroups[legIndex].add(read);
        }
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

    public List<SAMRecord> popCompletedReads(final ConsensusReads consensusReads)
    {
        List<SAMRecord> reads = Lists.newArrayList();

        for(int i = 0; i < mReadGroups.length; ++i)
        {
            if(mReadGroupComplete[i])
                continue;

            List<SAMRecord> readGroup = mReadGroups[i];

            if(readGroup == null || readGroup.size() < mFragmentCount)
                continue;

            if(mFragments.isEmpty()) // to avoid writing cached fragment reads twice
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

        if(!mFragments.isEmpty())
            mFragments.clear();

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

        for(List<SAMRecord> readGroup : mReadGroups)
        {
            if(readGroup != null && readGroup.size() == mFragmentCount)
            {
                return readGroup.stream().map(x -> x.getReadName()).collect(Collectors.toList());
            }
        }

        return null;
    }

    private enum ReadLegType
    {
        PRIMARY,
        MATE,
        PRIMARY_SUPPLEMENTARY,
        MATE_SUPPLEMENTARY;
    }

    private static ReadLegType getLegType(final SAMRecord read)
    {
        if(read.getReadPairedFlag() && !read.getFirstOfPairFlag())
        {
            return read.getSupplementaryAlignmentFlag() ? ReadLegType.MATE_SUPPLEMENTARY : ReadLegType.MATE;
        }
        else
        {
            return read.getSupplementaryAlignmentFlag() ? ReadLegType.PRIMARY_SUPPLEMENTARY : ReadLegType.PRIMARY;
        }
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

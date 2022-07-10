package com.hartwig.hmftools.svprep.reads;

import java.util.List;
import java.util.Set;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class ReadGroup
{
    private final List<ReadRecord> mReads;

    private GroupStatus mStatus;
    private boolean mSpansPartitions; // reads span a partition
    private Set<Integer> mJunctionPositions;

    public ReadGroup(final ReadRecord read)
    {
        mReads = Lists.newArrayListWithCapacity(2);
        mStatus = GroupStatus.UNSET;
        mSpansPartitions = false;
        mJunctionPositions = null;
        addRead(read);
    }

    public final String id() { return mReads.get(0).id(); }
    public List<ReadRecord> reads() { return mReads; }

    public boolean isComplete() { return mStatus == GroupStatus.COMPLETE; }
    public boolean isIncomplete() { return mStatus == GroupStatus.INCOMPLETE; }
    public boolean spansPartitions() { return mSpansPartitions; }
    public Set<Integer> junctionPositions() { return mJunctionPositions; }

    private enum GroupStatus
    {
        UNSET,
        INCOMPLETE,
        PARTIAL, // has all reads for initial partial, others are remote
        COMPLETE;
    }

    public void addRead(final ReadRecord read)
    {
        mReads.add(read);
    }

    public void addJunctionPosition(int position)
    {
        if(mJunctionPositions == null)
            mJunctionPositions = Sets.newHashSet();

        mJunctionPositions.add(position);
    }

    public String groupStatus() { return mStatus.toString(); }

    public void merge(final ReadGroup other)
    {
        other.reads().forEach(x -> addRead(x));
        setGroupStatus(null);
    }

    public boolean isSimpleComplete()
    {
        // no supplementaries and both reads received
        return mReads.size() == 2 && mReads.stream().allMatch(x -> !x.hasSuppAlignment() && !x.isSupplementaryAlignment());
    }

    public boolean allNoSupport() { return mReads.stream().allMatch(x -> x.readType() == ReadType.NO_SUPPORT); }

    public void setGroupStatus(final ChrBaseRegion region)
    {
        // use the provided partition region to determine if any reads in other partitions are missing
        boolean firstHasSupp = false;
        boolean secondHasSupp = false;

        for(ReadRecord read : mReads)
        {
            if(read.isFirstOfPair())
            {
                if(read.hasSuppAlignment() || read.isSupplementaryAlignment())
                    firstHasSupp = true;
            }
            else
            {
                if(read.hasSuppAlignment() || read.isSupplementaryAlignment())
                    secondHasSupp = true;
            }

            if(!mSpansPartitions && region != null)
            {
                if(read.hasSuppAlignment() && !supplementaryInRegion(read.supplementaryAlignment(), region))
                    mSpansPartitions = true;
                else if(!region.containsPosition(read.MateChromosome, read.MatePosStart))
                    mSpansPartitions = true;
            }
        }

        if(mSpansPartitions && region != null)
        {
            mStatus = GroupStatus.PARTIAL;
        }
        else
        {
            int requiredCount = 2 + (firstHasSupp ? 1 : 0) + (secondHasSupp ? 1 : 0);
            if(mReads.size() >= requiredCount)
                mStatus = GroupStatus.COMPLETE;
            else
                mStatus = GroupStatus.INCOMPLETE;
        }
    }

    private static boolean supplementaryInRegion(final SupplementaryReadData suppData, final ChrBaseRegion region)
    {
        return suppData != null && region.containsPosition(suppData.Chromosome, suppData.Position);
    }

    public int size() { return mReads.size(); }

    public String toString()
    {
        return String.format("%s reads(%d) state(%s) hasExternal(%s)", id(), mReads.size(), mStatus, mSpansPartitions);
    }
}

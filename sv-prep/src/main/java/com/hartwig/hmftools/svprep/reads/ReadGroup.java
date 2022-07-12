package com.hartwig.hmftools.svprep.reads;

import static com.hartwig.hmftools.svprep.CombinedReadGroups.formChromosomePartition;

import java.util.List;
import java.util.Set;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.CombinedReadGroups;

public class ReadGroup
{
    private final List<ReadRecord> mReads;

    private ReadGroupStatus mStatus;
    private boolean mSpansPartitions; // reads span a partition
    private Set<Integer> mJunctionPositions;

    public ReadGroup(final ReadRecord read)
    {
        mReads = Lists.newArrayListWithCapacity(2);
        mStatus = ReadGroupStatus.UNSET;
        mSpansPartitions = false;
        mJunctionPositions = null;
        addRead(read);
    }

    public final String id() { return mReads.get(0).id(); }
    public List<ReadRecord> reads() { return mReads; }

    public boolean isComplete() { return mStatus == ReadGroupStatus.COMPLETE; }
    public boolean isIncomplete() { return mStatus == ReadGroupStatus.INCOMPLETE; }
    public boolean spansPartitions() { return mSpansPartitions; }
    public Set<Integer> junctionPositions() { return mJunctionPositions; }

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

    public ReadGroupStatus groupStatus() { return mStatus; }

    public boolean isSimpleComplete()
    {
        // no supplementaries and both reads received
        return mReads.size() == 2 && mReads.stream().allMatch(x -> !x.hasSuppAlignment() && !x.isSupplementaryAlignment());
    }

    public boolean allNoSupport() { return mReads.stream().allMatch(x -> x.readType() == ReadType.NO_SUPPORT); }

    public ReadGroupState formGroupState(final ChrBaseRegion region, final String chrPartition, int partitionSize)
    {
        // use the provided partition region to determine if any reads in other partitions are missing
        boolean firstHasSupp = false;
        boolean secondHasSupp = false;
        String remoteChr = "";
        int remotePosition = 0;

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

            if(!mSpansPartitions)
            {
                if(read.hasSuppAlignment() && !supplementaryInRegion(read.supplementaryAlignment(), region))
                {
                    remoteChr = read.supplementaryAlignment().Chromosome;
                    remotePosition = read.supplementaryAlignment().Position;
                    mSpansPartitions = true;
                }
                else if(HumanChromosome.contains(read.MateChromosome) && !region.containsPosition(read.MateChromosome, read.MatePosStart))
                {
                    mSpansPartitions = true;
                    remoteChr = read.MateChromosome;
                    remotePosition = read.MatePosStart;
                }
            }
        }

        if(mSpansPartitions && region != null)
        {
            mStatus = ReadGroupStatus.PARTIAL;
        }
        else
        {
            int requiredCount = 2 + (firstHasSupp ? 1 : 0) + (secondHasSupp ? 1 : 0);
            if(mReads.size() >= requiredCount)
                mStatus = ReadGroupStatus.COMPLETE;
            else
                mStatus = ReadGroupStatus.INCOMPLETE;
        }

        String remoteChrPartition = remoteChr.isEmpty() ? "" : formChromosomePartition(remoteChr, remotePosition, partitionSize);

        return new ReadGroupState(
                id(), chrPartition, remoteChrPartition, remotePosition, mReads.size(), mStatus, firstHasSupp, secondHasSupp);
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

    /*
    public void merge(final ReadGroup other)
    {
        other.reads().forEach(x -> addRead(x));
        setGroupStatus(null);
    }
    */
}

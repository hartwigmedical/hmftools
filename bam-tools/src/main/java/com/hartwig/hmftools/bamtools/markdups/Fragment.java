package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.SUPPLEMENTARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNSET;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.formChromosomePartition;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.getFragmentCoordinates;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import htsjdk.samtools.SAMRecord;

public class Fragment
{
    private FragmentStatus mStatus;
    private final List<SAMRecord> mReads;
    private final FragmentCoordinates mCoordinates;
    private final boolean mUnpaired;
    private boolean mAllReadsPresent;
    private boolean mAllPrimaryReadsPresent;
    private List<String> mRemotePartitions; // partitions outside of the initial read's partition
    private boolean mReadsWritten;

    private double mAverageBaseQual;
    private int mDuplicateCount;

    public Fragment(final SAMRecord read)
    {
        mUnpaired = !read.getReadPairedFlag();

        if(!read.getSupplementaryAlignmentFlag())
        {
            mCoordinates = getFragmentCoordinates(read);
            mStatus = UNSET;
            mAllPrimaryReadsPresent = mUnpaired;
            mAllReadsPresent = mUnpaired && !read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE);
        }
        else
        {
            mCoordinates = null; // don't bother working it out
            mStatus = SUPPLEMENTARY;
            mAllPrimaryReadsPresent = false;
            mAllReadsPresent = false;
        }

        mReads = Lists.newArrayListWithCapacity(1); // most will only store the initial read
        mAverageBaseQual = 0;
        mDuplicateCount = 0;
        mRemotePartitions = null;
        mReads.add(read);
        mReadsWritten = false;
    }

    public final String id() { return mReads.get(0).getReadName(); }
    public List<SAMRecord> reads() { return mReads; }

    public FragmentStatus status() { return mStatus; }
    public void setStatus(final FragmentStatus status) { mStatus = status; }

    public boolean readsWritten() { return mReadsWritten; }
    public void setReadWritten() { mReadsWritten = true; }

    public boolean unpaired() { return mUnpaired; }
    public boolean allReadsPresent() { return mAllReadsPresent; }
    public boolean primaryReadsPresent() { return mAllPrimaryReadsPresent; }

    public FragmentCoordinates coordinates() { return mCoordinates; }
    public int initialPosition() { return mCoordinates.InitialPosition; }

    public double averageBaseQual() { return mAverageBaseQual; }
    public void setAverageBaseQual(double qual) { mAverageBaseQual = qual; }

    public int duplicateCount() { return mDuplicateCount; }
    public void setDuplicateCount(int count) { mDuplicateCount = count; }

    public void addRead(final SAMRecord read)
    {
        mReads.add(read);
        checkComplete();
    }

    public void merge(final Fragment other)
    {
        if(mReadsWritten)
        {
            BM_LOGGER.error("fragment({}) adding new fragment({}) when already written", this, other);
        }

        for(SAMRecord read : other.reads())
        {
            mReads.add(read);
        }

        if(other.hasRemotePartitions())
        {
            if(mRemotePartitions == null)
                mRemotePartitions = Lists.newArrayList();

            other.remotePartitions().stream().filter(x -> !mRemotePartitions.contains(x)).forEach(x -> mRemotePartitions.add(x));
        }

        checkComplete();
    }

    public boolean hasRemotePartitions() { return mRemotePartitions != null; }
    public List<String> remotePartitions() { return mRemotePartitions; }

    public void setRemotePartitions(final BaseRegion currentPartition)
    {
        mRemotePartitions = null;

        Set<String> chrPartitions = Sets.newHashSet();
        String chromosome = mReads.get(0).getContig();
        int partitionSize = currentPartition.baseLength();

        for(SAMRecord read : mReads)
        {
            /*
            if(read.getReadPairedFlag() && !read.getMateUnmappedFlag())
            {
                if(!read.getMateReferenceName().equals(chromosome) || !currentPartition.containsPosition(read.getMateAlignmentStart()))
                {
                    chrPartitions.add(formChromosomePartition(read.getMateReferenceName(), read.getMateAlignmentStart(), partitionSize));
                }
            }
            */

            if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
            {
                SupplementaryReadData suppData = SupplementaryReadData.from(read);

                if(!suppData.Chromosome.equals(chromosome) || !currentPartition.containsPosition(suppData.Position))
                {
                    chrPartitions.add(formChromosomePartition(suppData.Chromosome, suppData.Position, partitionSize));
                }
            }
        }

        if(!chrPartitions.isEmpty())
            mRemotePartitions = chrPartitions.stream().collect(Collectors.toList());
    }

    private void checkComplete()
    {
        if(mAllReadsPresent)
            return;

        int suppCount = 0;
        int nonSuppCount = 0;
        int expectedSuppCount = 0;
        int expectedNonSuppCount = 1;

        for(SAMRecord read : mReads)
        {
            if(read.getReadPairedFlag() && !read.getMateUnmappedFlag())
            {
                expectedNonSuppCount = 2;
            }

            if(read.getSupplementaryAlignmentFlag())
            {
                ++suppCount;
            }
            else
            {
                ++nonSuppCount;
            }

            if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
            {
                if(!read.getSupplementaryAlignmentFlag())
                {
                    ++expectedSuppCount;
                }
            }
        }

        mAllPrimaryReadsPresent = expectedNonSuppCount == nonSuppCount;
        mAllReadsPresent = (expectedNonSuppCount == nonSuppCount) && (expectedSuppCount == suppCount);
    }

    public String toString()
    {
        return String.format("reads(%d) status(%s) coords(%s) present(%s) id(%s) remotePartitions(%s)",
                mReads.size(), mStatus, mCoordinates.Key, mAllReadsPresent ? "all" : (mAllPrimaryReadsPresent ? "primary" : "incomplete"),
                id(), mRemotePartitions != null ? mRemotePartitions : "none");
    }
}

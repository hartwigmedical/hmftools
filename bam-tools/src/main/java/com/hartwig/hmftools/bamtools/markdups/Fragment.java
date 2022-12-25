package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.SUPPLEMENTARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNSET;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.formChromosomePartition;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.getUnclippedPosition;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.readToString;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.orientation;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import htsjdk.samtools.SAMRecord;

public class Fragment
{
    private FragmentStatus mStatus;
    private final List<SAMRecord> mReads;
    private final int[] mCoordinates; // negative/reverse strand positions are negated
    private final boolean mUnpaired;
    private boolean mAllReadsPresent;
    private boolean mAllPrimaryReadsPresent;
    private List<String> mRemotePartitions; // partitions outside of the initial read's partition
    private boolean mReadsWritten;

    public Fragment(final SAMRecord read)
    {
        mCoordinates = new int[SE_PAIR];
        mUnpaired = !read.getReadPairedFlag();

        if(!read.getSupplementaryAlignmentFlag())
        {
            int position = getUnclippedPosition(read);
            int strandPosition = orientation(read) == POS_ORIENT ? position : -position;
            mCoordinates[SE_START] = strandPosition;
            mStatus = UNSET;
            mAllPrimaryReadsPresent = mUnpaired;
            mAllReadsPresent = mUnpaired && !read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE);
        }
        else
        {
            mStatus = SUPPLEMENTARY;
            mAllPrimaryReadsPresent = false;
            mAllReadsPresent = false;
        }

        mReads = Lists.newArrayListWithCapacity(2);
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

    public int[] coordinates() { return mCoordinates; }
    public int initialPosition() { return mCoordinates[SE_START]; }

    public void addRead(final SAMRecord read)
    {
        mReads.add(read);

        if(mReadsWritten)
        {
            BM_LOGGER.error("fragment({}) adding new read({}) when already written", this, readToString(read));
        }

        checkComplete();

        if(!read.getSupplementaryAlignmentFlag())
        {
            int position = getUnclippedPosition(read);
            int strandPosition = orientation(read) == POS_ORIENT ? position : -position;

            if(abs(mCoordinates[SE_START]) <= position)
            {
                mCoordinates[SE_END] = strandPosition;
            }
            else
            {
                mCoordinates[SE_END] = mCoordinates[SE_START];
                mCoordinates[SE_START] = strandPosition;
            }
        }
    }

    public boolean hasRemotePartitions() { return mRemotePartitions != null; }
    public List<String> remotePartitions() { return mRemotePartitions; }

    public void setRemotePartitions(final BaseRegion currentPartition)
    {
        mRemotePartitions = null;

        // return true if links to other chromosomes
        List<String> chrPartitions = Lists.newArrayList();
        String chromosome = mReads.get(0).getContig();
        int partitionSize = currentPartition.baseLength();

        for(SAMRecord read : mReads)
        {
            if(read.getReadPairedFlag() && !read.getMateUnmappedFlag())
            {
                if(!read.getMateReferenceName().equals(chromosome) || !currentPartition.containsPosition(read.getMateAlignmentStart()))
                {
                    chrPartitions.add(formChromosomePartition(read.getMateReferenceName(), read.getMateAlignmentStart(), partitionSize));
                }
            }

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
            mRemotePartitions = chrPartitions;
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

    public int readCount() { return mReads.size(); }

    public String toString()
    {
        return String.format("reads(%d) status(%s) coords(%s:%d/%d) present(%s) id(%s) remotePartitions(%d)",
                mReads.size(), mStatus, mReads.get(0).getContig(), mCoordinates[SE_START], mCoordinates[SE_END],
                mAllReadsPresent ? "all" : (mAllPrimaryReadsPresent ? "primary" : "incomplete"), id(),
                mRemotePartitions != null ? mRemotePartitions.size() : 0);
    }
}

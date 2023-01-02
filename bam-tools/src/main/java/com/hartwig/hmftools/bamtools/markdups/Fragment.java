package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.markdups.FragmentCoordinates.NO_COORDS;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.SUPPLEMENTARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNSET;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.formChromosomePartition;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.getFragmentCoordinates;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import htsjdk.samtools.SAMRecord;

public class Fragment
{
    private FragmentStatus mStatus;
    private final boolean mUnpaired;
    private boolean mAllReadsPresent;
    private boolean mAllPrimaryReadsPresent;
    private final List<SAMRecord> mReads; // consider making an array of 4 (or less for BNDs)
    private FragmentCoordinates mCoordinates;

    // private String mBasePartition; // partition of the lowest non-supplementary read

    @Deprecated
    private List<String> mRemotePartitions; // partitions outside of the initial read's partition. currently only used for supplementaries

    private boolean mReadsWritten;

    private double mAverageBaseQual;
    private int mDuplicateCount;

    public Fragment(final SAMRecord read)
    {
        mUnpaired = !read.getReadPairedFlag();

        mReads = Lists.newArrayListWithCapacity(2);
        mReads.add(read);

        // mBasePartition = null;

        if(!read.getSupplementaryAlignmentFlag())
        {
            // mCoordinates = getFragmentCoordinates(read, false);
            mStatus = UNSET;

            if(mUnpaired)
            {
                mAllPrimaryReadsPresent = true;
                mAllReadsPresent = !read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE);
            }
            else
            {
                mAllPrimaryReadsPresent = false;
                mAllReadsPresent = false;
            }

            mCoordinates = getFragmentCoordinates(mReads);
        }
        else
        {
            mCoordinates = NO_COORDS; // don't bother working it out
            mStatus = SUPPLEMENTARY;
            mAllPrimaryReadsPresent = false;
            mAllReadsPresent = false;
        }

        mAverageBaseQual = 0;
        mDuplicateCount = 0;
        mRemotePartitions = null;
        mReadsWritten = false;
    }

    public final String id() { return mReads.get(0).getReadName(); }
    public List<SAMRecord> reads() { return mReads; }

    public FragmentStatus status() { return mStatus; }
    public void setStatus(final FragmentStatus status) { mStatus = status; }

    public boolean unpaired() { return mUnpaired; }
    public boolean allReadsPresent() { return mAllReadsPresent; }
    public boolean primaryReadsPresent() { return mAllPrimaryReadsPresent; }

    public FragmentCoordinates coordinates() { return mCoordinates; }
    public int initialPosition() { return mCoordinates.InitialPosition; }

    public double averageBaseQual() { return mAverageBaseQual; }
    public void setAverageBaseQual(double qual) { mAverageBaseQual = qual; }

    public int duplicateCount() { return mDuplicateCount; }
    public void setDuplicateCount(int count) { mDuplicateCount = count; }

    public boolean readsWritten() { return mReadsWritten; }
    public void setReadWritten() { mReadsWritten = true; }

    public void addRead(final SAMRecord read)
    {
        mReads.add(read);

        if(!read.getSupplementaryAlignmentFlag())
        {
            mAllPrimaryReadsPresent = true;
            mCoordinates = getFragmentCoordinates(mReads);
        }

        checkComplete();
    }

    // public String basePartition() { return mBasePartition; }
    // public void setBasePartition(final String partition) { mBasePartition = partition; }

    public boolean hasRemotePartitions() { return mRemotePartitions != null; }
    public List<String> remotePartitions() { return mRemotePartitions; }

    public static String getBasePartition(final SAMRecord read, final int partitionSize)
    {
        if(read.getSupplementaryAlignmentFlag())
        {
            SupplementaryReadData suppData = SupplementaryReadData.from(read);
            if(suppData != null)
                return formChromosomePartition(suppData.Chromosome, suppData.Position, partitionSize);
        }

        if(!read.getReadPairedFlag())
            return formChromosomePartition(read.getReferenceName(), read.getAlignmentStart(), partitionSize);

        // take the lower of the read and its mate
        boolean readLowerPos;
        if(read.getReferenceIndex() == read.getMateReferenceIndex())
        {
            readLowerPos = read.getAlignmentStart() < read.getMateAlignmentStart();
        }
        else
        {
            readLowerPos = read.getReferenceIndex() < read.getMateReferenceIndex();
        }

        return readLowerPos ?
                formChromosomePartition(read.getReferenceName(), read.getAlignmentStart(), partitionSize)
                : formChromosomePartition(read.getMateReferenceName(), read.getMateAlignmentStart(), partitionSize);
    }

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
                chrPartitions.add(formChromosomePartition(suppData.Chromosome, suppData.Position, partitionSize));
            }
        }

        if(!chrPartitions.isEmpty())
            mRemotePartitions = chrPartitions.stream().collect(Collectors.toList());
    }

    public int expectedReadCount()
    {
        int expectedSuppCount = 0;
        int expectedNonSuppCount = 1;

        for(SAMRecord read : mReads)
        {
            if(read.getReadPairedFlag() && !read.getMateUnmappedFlag())
            {
                expectedNonSuppCount = 2;
            }

            if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE) && !read.getSupplementaryAlignmentFlag())
            {
                ++expectedSuppCount;
            }
        }

        return expectedNonSuppCount + expectedSuppCount;
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
        return String.format("id(%s) reads(%d) status(%s) coords(%s) present(%s) remotePartitions(%s)",
                id(), mReads.size(), mStatus, mCoordinates.Key,
                mAllReadsPresent ? "all" : (mAllPrimaryReadsPresent ? "primary" : "incomplete"),
                mRemotePartitions != null ? mRemotePartitions : "none");
    }
}

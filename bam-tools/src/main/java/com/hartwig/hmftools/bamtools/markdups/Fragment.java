package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.NONE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.SUPPLEMENTARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNCLEAR;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNSET;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.formChromosomePartition;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.getUnclippedPosition;
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

    public boolean readWritten() { return mReadsWritten; }
    public void setReadWritten() { mReadsWritten = true; }

    public boolean unpaired() { return mUnpaired; }
    public boolean allReadsPresent() { return mAllReadsPresent; }
    public boolean primaryReadsPresent() { return mAllPrimaryReadsPresent; }

    public int[] coordinates() { return mCoordinates; }
    public int initialPosition() { return mCoordinates[SE_START]; }

    public void addRead(final SAMRecord read)
    {
        mReads.add(read);

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

    public static FragmentStatus calcFragmentStatus(final Fragment first, final Fragment second)
    {
        if(first.unpaired() != second.unpaired())
            return NONE;

        if(first.primaryReadsPresent() == second.primaryReadsPresent())
        {
            if(first.unpaired())
            {
                return first.initialPosition() == second.initialPosition() ? DUPLICATE : NONE;
            }
            else
            {
                return first.coordinates()[SE_START] == second.coordinates()[SE_START]
                    && first.coordinates()[SE_END] == second.coordinates()[SE_END] ? DUPLICATE : NONE;
            }
        }
        else
        {
            if(first.initialPosition() != second.initialPosition())
                return NONE;

            // mate start positions must be within close proximity
            SAMRecord firstRead = first.reads().get(0);
            SAMRecord secondRead = second.reads().get(0);

            if(!firstRead.getMateReferenceName().equals(secondRead.getMateReferenceName()))
                return NONE;

            return abs(firstRead.getMateAlignmentStart() - secondRead.getMateAlignmentStart()) < firstRead.getReadLength()
                    ? UNCLEAR : NONE;
        }
    }

    public boolean hasRemotePartitions() { return mRemotePartitions != null; }
    public List<String> remotePartitions() { return mRemotePartitions; }

    public boolean setRemotePartitions(final BaseRegion currentPartition)
    {
        // return true if links to other chromosomes
        List<String> chrPartitions = Lists.newArrayList();
        String chromosome = mReads.get(0).getContig();
        int partitionSize = currentPartition.baseLength();

        boolean hasRemoteChromosomes = false;

        for(SAMRecord read : mReads)
        {
            if(read.getReadPairedFlag() && !read.getMateUnmappedFlag())
            {
                if(!read.getMateReferenceName().equals(chromosome))
                {
                    hasRemoteChromosomes = true;

                    if(!currentPartition.containsPosition(read.getMateAlignmentStart()))
                    {
                        chrPartitions.add(formChromosomePartition(read.getMateReferenceName(), read.getMateAlignmentStart(), partitionSize));
                    }
                }
            }

            if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
            {
                SupplementaryReadData suppData = SupplementaryReadData.from(read);

                if(!suppData.Chromosome.equals(chromosome))
                {
                    hasRemoteChromosomes = true;

                    if(!currentPartition.containsPosition(suppData.Position))
                    {
                        chrPartitions.add(formChromosomePartition(suppData.Chromosome, suppData.Position, partitionSize));
                    }
                }
            }
        }

        if(!chrPartitions.isEmpty())
            mRemotePartitions = chrPartitions;

        return hasRemoteChromosomes;
    }

    private void checkComplete()
    {
        if(mAllReadsPresent)
            return;

        int suppCount = 0;
        int nonSuppCount = 0;
        int expectedSuppCount = 0;
        int expectedNonSuppCount = 0;

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

        mAllReadsPresent = (expectedNonSuppCount == nonSuppCount) && (expectedSuppCount == suppCount);
    }

    public int size() { return mReads.size(); }

    public String toString()
    {
        return String.format("reads(%d) status(%s) coords(%d - %d) initReadStart(%s:%d) id(%s)",
                mReads.size(), mStatus, mReads.get(0).getContig(), mReads.get(0).getAlignmentStart(), id());
    }


    /*

            BaseRegion currentRange = null;
        BaseRegion expectedRange = null;

        for(SAMRecord read : mReads)
        {
            if(read.getReadPairedFlag() && !read.getMateUnmappedFlag())
            {
                if(read.getMateReferenceName().equals(chromosome))
                {
                    expectedRange.setStart(min(expectedRange.start(), read.getMateAlignmentStart()));
                    expectedRange.setEnd(max(expectedRange.end(), read.getMateAlignmentStart()));
                }
                else
                {
                    chrPartitions.add(formChromosomePartition(read.getMateReferenceName(), read.getMateAlignmentStart(), partitionSize));
                }
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

                SupplementaryReadData suppData = SupplementaryReadData.from(read);

                if((suppData.Chromosome.equals(chromosome)) && currentPartition.containsPosition(suppData.Position))
                {
                    expectedRange.setStart(min(expectedRange.start(), suppData.Position));
                    expectedRange.setEnd(max(expectedRange.end(), suppData.Position));
                }
                else
                {
                    chrPartitions.add(formChromosomePartition(suppData.Chromosome, suppData.Position, partitionSize));
                }
            }
        }

        boolean isComplete = (expectedNonSuppCount == nonSuppCount) && (expectedSuppCount == suppCount);

     */

}

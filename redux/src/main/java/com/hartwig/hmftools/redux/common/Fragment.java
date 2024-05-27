package com.hartwig.hmftools.redux.common;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAP_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getFivePrimeUnclippedPosition;
import static com.hartwig.hmftools.redux.common.Constants.DEFAULT_POS_BUFFER_SIZE;
import static com.hartwig.hmftools.redux.common.FragmentCoordinates.NO_COORDS;
import static com.hartwig.hmftools.redux.common.FragmentStatus.SUPPLEMENTARY;
import static com.hartwig.hmftools.redux.common.FragmentStatus.UNSET;
import static com.hartwig.hmftools.redux.common.FragmentUtils.formChromosomePartition;
import static com.hartwig.hmftools.redux.common.FragmentUtils.getFragmentCoordinates;
import static com.hartwig.hmftools.redux.common.ReadUnmapper.parseUnmappedCoords;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

public class Fragment
{
    private FragmentStatus mStatus;
    private final boolean mUnpaired;
    private final boolean mHasLocalMate;
    private boolean mAllReadsPresent;
    private boolean mAllPrimaryReadsPresent;
    private final List<SAMRecord> mReads; // consider making an array of 4 (or less for BNDs)
    private FragmentCoordinates mCoordinates;

    private boolean mReadsWritten;

    // duplicate read info
    private double mAverageBaseQual;
    private String mCandidateDupKey;
    private String mUmi;

    public Fragment(final SAMRecord read)
    {
        mUnpaired = !read.getReadPairedFlag();

        mReads = Lists.newArrayListWithCapacity(2);
        mReads.add(read);

        if(!read.getSupplementaryAlignmentFlag())
        {
            mStatus = UNSET;

            if(mUnpaired)
            {
                mAllPrimaryReadsPresent = true;
                mAllReadsPresent = !read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE);
                mHasLocalMate = false;
            }
            else
            {
                mAllPrimaryReadsPresent = false;
                mAllReadsPresent = false;

                if(read.getMateUnmappedFlag()) // unmapped reads come through the same slice
                {
                    String mateCoordsStr = read.getStringAttribute(UNMAP_ATTRIBUTE);
                    if(mateCoordsStr != null)
                    {
                        String[] mateCoords = parseUnmappedCoords(mateCoordsStr);
                        String mateChr = mateCoords[0];
                        int matePosition = Integer.parseInt(mateCoords[1]);

                        mHasLocalMate = mateChr.equals(read.getReferenceName())
                                && abs(matePosition - read.getAlignmentStart()) < DEFAULT_POS_BUFFER_SIZE;
                    }
                    else
                    {
                        mHasLocalMate = true;
                    }
                }
                else
                {
                    mHasLocalMate = read.getMateReferenceName().equals(read.getReferenceName())
                            && abs(read.getMateAlignmentStart() - read.getAlignmentStart()) < DEFAULT_POS_BUFFER_SIZE;
                }
            }
        }
        else
        {
            mStatus = SUPPLEMENTARY;
            mAllPrimaryReadsPresent = false;
            mAllReadsPresent = false;
            mHasLocalMate = false;
        }

        mCoordinates = NO_COORDS; // unset for non primary reads

        mAverageBaseQual = 0;
        mReadsWritten = false;
        mCandidateDupKey = null;
        mUmi = null;
    }

    public final String id() { return mReads.get(0).getReadName(); }
    public List<SAMRecord> reads() { return mReads; }

    public FragmentStatus status() { return mStatus; }
    public void setStatus(final FragmentStatus status) { mStatus = status; }

    public boolean unpaired() { return mUnpaired; }
    public boolean allReadsPresent() { return mAllReadsPresent; }
    public boolean primaryReadsPresent() { return mAllPrimaryReadsPresent; }
    public boolean hasLocalMate() { return mHasLocalMate; }

    public FragmentCoordinates coordinates() { return mCoordinates; }
    public void setCoordinates(final FragmentCoordinates coordinates) { mCoordinates = coordinates; }
    public int initialPosition() { return mCoordinates.InitialPosition; }

    public void intialiseCoordinates(boolean useMateCigar) { mCoordinates = getFragmentCoordinates(mReads, useMateCigar); }

    public double averageBaseQual() { return mAverageBaseQual; }
    public void setAverageBaseQual(double qual) { mAverageBaseQual = qual; }

    public String candidateDupKey() { return mCandidateDupKey; }
    public void setCandidateDupKey(final String key) { mCandidateDupKey = key; }

    public String umi() { return mUmi; }
    public void setUmi(final String umi) { mUmi = umi; }

    public boolean readsWritten() { return mReadsWritten; }
    public void setReadWritten() { mReadsWritten = true; }

    public void addRead(final SAMRecord read)
    {
        mReads.add(read);

        if(!read.getSupplementaryAlignmentFlag())
        {
            mAllPrimaryReadsPresent = true;

            if(mCoordinates.Incomplete)
                mCoordinates = getFragmentCoordinates(mReads, false);
        }

        checkComplete();
    }

    public static String getBasePartition(final SAMRecord read, final int partitionSize)
    {
        if(read.getSupplementaryAlignmentFlag())
        {
            // get the lower of their mate or their supplementary read's location
            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read);

            boolean hasSuppData = suppData != null && HumanChromosome.contains(suppData.Chromosome);

            boolean hasMate = read.getReadPairedFlag() && !read.getMateUnmappedFlag() && HumanChromosome.contains(read.getMateReferenceName());

            if(!hasSuppData && !hasMate)
                return null;

            boolean useMate;

            if(hasSuppData && hasMate)
            {
                if(suppData.Chromosome.equals(read.getMateReferenceName()))
                    useMate = read.getMateAlignmentStart() <= suppData.Position;
                else
                    useMate = HumanChromosome.lowerChromosome(read.getMateReferenceName(), suppData.Chromosome);
            }
            else
            {
                useMate = hasMate;
            }

            return useMate ?
                    formChromosomePartition(read.getMateReferenceName(), read.getMateAlignmentStart(), partitionSize)
                    : formChromosomePartition(suppData.Chromosome, suppData.Position, partitionSize);
        }

        if(!read.getReadPairedFlag())
            return formChromosomePartition(read.getReferenceName(), read.getAlignmentStart(), partitionSize);

        if(!HumanChromosome.contains(read.getMateReferenceName()))
            return null;

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

                if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
                {
                    ++expectedSuppCount;
                }
            }
        }

        mAllPrimaryReadsPresent = expectedNonSuppCount == nonSuppCount;
        mAllReadsPresent = (expectedNonSuppCount == nonSuppCount) && (expectedSuppCount == suppCount);
    }

    public int readCount() { return mReads.size(); }

    public boolean isPreciseInversion()
    {
        SAMRecord first = mReads.stream().filter(x -> !x.getSupplementaryAlignmentFlag() && x.getFirstOfPairFlag()).findFirst().orElse(null);
        SAMRecord second = mReads.stream().filter(x -> !x.getSupplementaryAlignmentFlag() && x.getSecondOfPairFlag()).findFirst().orElse(null);

        if(first == null || second == null)
            return false;

        return first.getReadNegativeStrandFlag() == second.getReadNegativeStrandFlag()
                && !first.getReadUnmappedFlag() && !second.getReadUnmappedFlag()
                && getFivePrimeUnclippedPosition(first) == getFivePrimeUnclippedPosition(second);
    }

    public String toString()
    {
        return String.format("id(%s) reads(%d) status(%s) coords(%s) present(%s)", //  mate(%s:%d)
                id(), mReads.size(), mStatus, mCoordinates.keyOriented(),
                mAllReadsPresent ? "all" : (mAllPrimaryReadsPresent ? "primary" : "incomplete"));
    }
}

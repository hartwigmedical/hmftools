package com.hartwig.hmftools.sieve.annotate;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateNegativeStrand;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateUnmapped;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

// TODO(m_cooper): Duplication.
public class AnnotatedBedRecord
{
    public static final String TSV_HEADER =
            "Chromosome\tPosStart\tPosEnd\tSampleCount\tDepthMin\tDepthMax\tPrimaryReadCount\tPrimarySoftClippedCount\tSupplementaryCount\tPrimaryImproperPairCount";

    private final String mChromosome;
    private final int mPosStart;
    private final int mPosEnd;
    private final long mSampleCount;
    private final long mDepthMin;
    private final long mDepthMax;

    private long mPrimaryReadCount;
    private long mPrimarySoftClippedCount;
    private long mSupplementaryCount;
    private long mPrimaryImproperPairCount;

    public AnnotatedBedRecord(@NotNull final String chromosome, final int posStart, final int posEnd, final long sampleCount,
            final long depthMin, final long depthMax)
    {
        mChromosome = stripChrPrefix(chromosome);
        mPosStart = posStart;
        mPosEnd = posEnd;
        mSampleCount = sampleCount;
        mDepthMin = depthMin;
        mDepthMax = depthMax;

        mPrimaryReadCount = 0;
        mPrimarySoftClippedCount = 0;
        mSupplementaryCount = 0;
        mPrimaryImproperPairCount = 0;
    }

    // TODO(m_cooper): Use the original implementation.
    private static boolean isNotProperReadPair(@NotNull final SAMRecord read)
    {
        if(read.getReadUnmappedFlag())
        {
            return false;
        }

        // or a fragment length outside the observed distribution
        // TODO(m_cooper): Make this configurable.
        if(abs(read.getInferredInsertSize()) > 1000)
        {
            return true;
        }

        // an unmapped mate
        if(mateUnmapped(read))
        {
            return true;
        }

        if(read.getReadPairedFlag())
        {
            // inter-chromosomal
            if(!read.getReferenceName().equals(read.getMateReferenceName()))
            {
                return true;
            }

            // inversion
            return read.getReadNegativeStrandFlag() == mateNegativeStrand(read);
        }

        return false;
    }

    public void matchedRead(@NotNull final SAMRecord read)
    {
        if(read.getSupplementaryAlignmentFlag())
        {
            mSupplementaryCount++;
            return;
        }

        mPrimaryReadCount++;

        final Cigar cigar = read.getCigar();
        final var it = cigar.getCigarElements().iterator();
        while(it.hasNext())
        {
            final CigarElement el = it.next();
            final CigarOperator op = el.getOperator();
            if(op == CigarOperator.SOFT_CLIP)
            {
                mPrimarySoftClippedCount++;
                break;
            }
        }

        if(isNotProperReadPair(read))
        {
            // TODO(m_cooper): Double counting?
            mPrimaryImproperPairCount++;
        }
    }

    public void resetCounts()
    {
        mPrimaryReadCount = 0;
        mPrimarySoftClippedCount = 0;
        mSupplementaryCount = 0;
        mPrimaryImproperPairCount = 0;
    }

    @NotNull
    public String getChromosome()
    {
        return mChromosome;
    }

    public int getPosStart()
    {
        return mPosStart;
    }

    public int getPosEnd()
    {
        return mPosEnd;
    }

    public long getSampleCount()
    {
        return mSampleCount;
    }

    public long getDepthMin()
    {
        return mDepthMin;
    }

    public long getDepthMax()
    {
        return mDepthMax;
    }

    public long getReadCount()
    {
        return mPrimaryReadCount;
    }

    public long getSoftClippedCount()
    {
        return mPrimarySoftClippedCount;
    }

    public long getSupplementaryCount()
    {
        return mSupplementaryCount;
    }

    public long getImproperPairCount()
    {
        return mPrimaryImproperPairCount;
    }

    @Override
    public String toString()
    {
        final String sb = mChromosome + '\t' + mPosStart + '\t' + mPosEnd + '\t' + mSampleCount + '\t' + mDepthMin + '\t' + mDepthMax + '\t'
                + mPrimaryReadCount + '\t' + mPrimarySoftClippedCount + '\t' + mSupplementaryCount + '\t' + mPrimaryImproperPairCount;

        return sb;
    }
}

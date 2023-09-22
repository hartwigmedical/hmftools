package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;

import java.util.Objects;

import org.jetbrains.annotations.NotNull;

public class BlacklistRegion
{
    public static final String CSV_HEADER = "Chromosome,PosStart,PosEnd,SampleCount,DepthMin,DepthMax";

    private final String mChromosome;
    private final int mPosStart;
    private final int mPosEnd;
    private final int mSampleCount;
    private final int mDepthMin;
    private final int mDepthMax;

    public BlacklistRegion(@NotNull final String chromosome, final int posStart, final int posEnd, final int sampleCount,
            final int depthMin, final int depthMax)
    {
        mChromosome = stripChrPrefix(chromosome);
        mPosStart = posStart;
        mPosEnd = posEnd;
        mSampleCount = sampleCount;
        mDepthMin = depthMin;
        mDepthMax = depthMax;
    }

    public String getCSVFragment()
    {
        return mChromosome
                + ','
                + mPosStart
                + ','
                + mPosEnd
                + ','
                + mSampleCount
                + ','
                + mDepthMin
                + ','
                + mDepthMax;
    }

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

    public int getSampleCount()
    {
        return mSampleCount;
    }

    public int getDepthMin()
    {
        return mDepthMin;
    }

    public int getDepthMax()
    {
        return mDepthMax;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
        {
            return true;
        }
        if(!(o instanceof BlacklistRegion))
        {
            return false;
        }
        final BlacklistRegion that = (BlacklistRegion) o;
        return mPosStart == that.mPosStart && mPosEnd == that.mPosEnd && mSampleCount == that.mSampleCount && mDepthMin == that.mDepthMin
                && mDepthMax == that.mDepthMax && mChromosome.equals(that.mChromosome);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mChromosome, mPosStart, mPosEnd, mSampleCount, mDepthMin, mDepthMax);
    }

    @Override
    public String toString()
    {
        return "BlacklistRegion{" +
                "Chromosome='" + mChromosome + '\'' +
                ", PosStart=" + mPosStart +
                ", PosEnd=" + mPosEnd +
                ", SampleCount=" + mSampleCount +
                ", DepthMin=" + mDepthMin +
                ", DepthMax=" + mDepthMax +
                '}';
    }
}

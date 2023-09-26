package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;

public class BlacklistRegion
{
    public static final String TSV_HEADER = "Chromosome\tPosStart\tPosEnd\tSampleCount\tDepthMin\tDepthMax";

    private final String mChromosome;
    private final int mPosStart;
    private final int mPosEnd;
    private final int mSampleCount;
    private final int mDepthMin;
    private final int mDepthMax;

    public BlacklistRegion(final String chromosome, final int posStart, final int posEnd, final int sampleCount,
            final int depthMin, final int depthMax)
    {
        mChromosome = stripChrPrefix(chromosome);
        mPosStart = posStart;
        mPosEnd = posEnd;
        mSampleCount = sampleCount;
        mDepthMin = depthMin;
        mDepthMax = depthMax;
    }

    public String getTSVFragment()
    {
        return mChromosome
                + '\t'
                + mPosStart
                + '\t'
                + mPosEnd
                + '\t'
                + mSampleCount
                + '\t'
                + mDepthMin
                + '\t'
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
}

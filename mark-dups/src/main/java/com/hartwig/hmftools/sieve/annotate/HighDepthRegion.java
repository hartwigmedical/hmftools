package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;

public class HighDepthRegion
{
    public static final String TSV_HEADER = "Chromosome\tPosStart\tPosEnd\tBaseDepth";

    private final String mChromosome;
    private final int mPosStart;
    private final int mPosEnd;
    private final int mBaseDepth;

    public HighDepthRegion(final String chromosome, final int posStart, final int posEnd, final int baseDepth)
    {
        mChromosome = stripChrPrefix(chromosome);
        mPosStart = posStart;
        mPosEnd = posEnd;
        mBaseDepth = baseDepth;
    }

    public String getTSVFragment()
    {
        return mChromosome
                + '\t'
                + mPosStart
                + '\t'
                + mPosEnd
                + '\t'
                + mBaseDepth;
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
}

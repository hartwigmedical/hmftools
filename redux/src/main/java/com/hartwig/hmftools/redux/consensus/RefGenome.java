package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.redux.consensus.BaseQualPair.NO_BASE;

import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class RefGenome
{
    private final RefGenomeInterface mRefGenome;
    private final Map<String,Integer> mChromosomeLengths;

    public RefGenome(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
        mChromosomeLengths = mRefGenome.chromosomeLengths();
    }

    public byte getRefBase(final String chromosome, int position)
    {
        int chromosomeLength = mChromosomeLengths.getOrDefault(chromosome, 0);

        if(position > chromosomeLength)
            return NO_BASE;

        return mRefGenome.getBases(chromosome, position, position)[0];
    }

    public byte[] getRefBases(final String chromosome, int posStart, int posEnd)
    {
        int chromosomeLength = mChromosomeLengths.getOrDefault(chromosome, 0);

        if(posEnd > chromosomeLength)
            return null;

        return mRefGenome.getBases(chromosome, posStart, posEnd);
    }

    public int getChromosomeLength(final String chromosome) { return mChromosomeLengths.get(chromosome); }
}

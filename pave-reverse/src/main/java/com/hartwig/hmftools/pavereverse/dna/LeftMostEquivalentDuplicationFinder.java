package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class LeftMostEquivalentDuplicationFinder extends LeftMostEquivalentChangeFinder
{
    private final String mRef;

    public LeftMostEquivalentDuplicationFinder(RefGenomeInterface genome, String chromosome, int start, int end)
    {
        super(genome, chromosome, start, end - start);
        mRef = getBases(mStart);
    }

    @Override
    boolean changeAtPositionIsEquivalent(final int position)
    {
        String atPosition = getBases(position);
        return atPosition.equals(mRef);
    }

    private String getBases(int position)
    {
        return mGenome.getBaseString(mChromosome, position, position + mLength);
    }
}

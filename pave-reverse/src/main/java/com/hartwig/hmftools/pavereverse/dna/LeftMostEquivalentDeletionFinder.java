package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class LeftMostEquivalentDeletionFinder extends LeftMostEquivalentChangeFinder
{
    private final String mRef;

    public LeftMostEquivalentDeletionFinder(RefGenomeInterface genome, String chromosome, int start, int end)
    {
        super(genome, chromosome, start, end - start);
        mRef = getBases(mStart);
    }

    @Override
    boolean changeAtPositionIsEquivalent(final int position)
    {
        // S = vwxyz
        // S_del_y = S_del_w <=> vwxz = vxyz <=> vwx = vxy <=> wx = xy
        String fromCurrentPositionWithDeletionAtOriginalLocation = mGenome.getBaseString(mChromosome, position, mStart - 1);
        String withDeletionAtCurrentPosition = mGenome.getBaseString(mChromosome, position + mLength + 1,mStart + mLength);
        return fromCurrentPositionWithDeletionAtOriginalLocation.equals(withDeletionAtCurrentPosition);
    }

    private String getBases(int position)
    {
        return mGenome.getBaseString(mChromosome, position, position + mLength);
    }
}

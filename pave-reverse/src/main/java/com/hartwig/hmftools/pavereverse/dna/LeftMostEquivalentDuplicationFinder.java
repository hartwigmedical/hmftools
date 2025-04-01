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
        // S = vwxyz
        // S_dup(w) = S_dup(y) <=> vwwxyz = vwxyyz <=> vwwxy = vwxyy <=> vwwx = vwxy <=> wwx = wxy <=> wx = xy

        String basesToDuplicateAtCurrentPosition = getBases(position);
        String basesFromCurrentPositionToOriginalPosition = mGenome.getBaseString(mChromosome, position, mStart - 1);
        String withDuplicationAtOriginalPosition = basesFromCurrentPositionToOriginalPosition + mRef;
        String withDuplicationAtCurrentPosition = basesToDuplicateAtCurrentPosition + basesFromCurrentPositionToOriginalPosition;
        return withDuplicationAtOriginalPosition.equals(withDuplicationAtCurrentPosition);
    }

    private String getBases(int position)
    {
        return mGenome.getBaseString(mChromosome, position, position + mLength);
    }
}

package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public abstract class LeftMostEquivalentChangeFinder
{
    final RefGenomeInterface mGenome;
    final String mChromosome;
    final int mStart;
    final int mLength;

    public LeftMostEquivalentChangeFinder(RefGenomeInterface genome, String chromosome, int start, int length)
    {
        mGenome = genome;
        mChromosome = chromosome;
        mStart = start;
        mLength = length;
    }

    abstract boolean changeAtPositionIsEquivalent(int position);

    public int findLeftMostEquivalentPosition()
    {
        int position = mStart;
        int result = mStart;
        int distanceFromPreviousMatch = 0;
        while(position > 0)
        {
            if(changeAtPositionIsEquivalent(position))
            {
                result = position;
                distanceFromPreviousMatch = 0;
            }
            else
            {
                if(distanceFromPreviousMatch >= mLength)
                {
                    return result;
                }
                distanceFromPreviousMatch++;
            }
            position--;
        }
        return result;
    }
}

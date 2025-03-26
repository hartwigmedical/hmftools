package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class LeftMostEquivalentInsertionFinder extends LeftMostEquivalentChangeFinder
{
    private final String mToInsert;

    public LeftMostEquivalentInsertionFinder(RefGenomeInterface genome, String chromosome, int position, String bases)
    {
        super(genome, chromosome, position, bases.length());
        mToInsert = bases;
    }

    @Override
    boolean changeAtPositionIsEquivalent(final int position)
    {
        String basesBetweenCurrentPositionAndOriginalPosition = getBases(position, mStart);
        String baseAtCurrentPosition = basesBetweenCurrentPositionAndOriginalPosition.substring(0, 1);
        String remainder = basesBetweenCurrentPositionAndOriginalPosition.substring(1);
        String withInsertionImmediatelyAfterCurrentPosition = baseAtCurrentPosition + mToInsert + remainder;
        String withInsertionAtOriginalPosition =  basesBetweenCurrentPositionAndOriginalPosition + mToInsert;
        return withInsertionImmediatelyAfterCurrentPosition.equals(withInsertionAtOriginalPosition);
    }

    private String getBases(int start, int stop)
    {
        return mGenome.getBaseString(mChromosome, start, stop);
    }
}

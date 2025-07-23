package com.hartwig.hmftools.pavereverse.dna;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.pavereverse.FixedStringGenome;

import org.junit.Test;

public class LeftMostEquivalentInsertionFinderTest
{
    private RefGenomeInterface genome;

    @Test
    public void alreadyLeftMost()
    {
        genome = new FixedStringGenome("AACCTTGG");
        LeftMostEquivalentInsertionFinder finder = finder( 4, "G");
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder( 3, "G");
        assertEquals(3, finder.findLeftMostEquivalentPosition());
    }

    @Test
    public void repeatedBaseInsertion()
    {
        genome = new FixedStringGenome("CCTGAAAAAAGG");
        LeftMostEquivalentInsertionFinder finder = finder( 4, "AAAAAAA");
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder( 5, "AAAAAAA");
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder( 6, "AAAAAAAA");
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(7, "AAAAAAA");
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder( 10, "AAAAAAA");
        assertEquals(4, finder.findLeftMostEquivalentPosition());
    }

    @Test
    public void length2RepeatPatternInsertion()
    {
        genome = new FixedStringGenome("ACCGCGCGTT");
        LeftMostEquivalentInsertionFinder finder = finder(3, "CGCGCGCG");
        assertEquals(3, finder.findLeftMostEquivalentPosition());

        finder = finder(4, "CGCGCGCG");
        assertEquals(2, finder.findLeftMostEquivalentPosition());

        finder = finder(6, "CGCGCGCG");
        assertEquals(2, finder.findLeftMostEquivalentPosition());

        finder = finder(8, "CGCGCGCG");
        assertEquals(2, finder.findLeftMostEquivalentPosition());
    }

    private LeftMostEquivalentInsertionFinder finder(int position, String bases)
    {
        return new LeftMostEquivalentInsertionFinder(genome, "1", position, bases);
    }
}

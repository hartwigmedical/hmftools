package com.hartwig.hmftools.pavereverse.dna;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.pavereverse.FixedStringGenome;

import org.junit.Test;

public class LeftMostEquivalentDeletionFinderTest
{
    private RefGenomeInterface genome;
    
    @Test
    public void singleBaseDeletion()
    {
        genome = new FixedStringGenome("CCTGAAAAAAGG");
        LeftMostEquivalentDeletionFinder finder = finder( 1, 1);
        assertEquals(1, finder.findLeftMostEquivalentPosition());

        finder = finder( 2, 2);
        assertEquals(1, finder.findLeftMostEquivalentPosition());

        finder = finder(3, 3);
        assertEquals(3, finder.findLeftMostEquivalentPosition());

        finder = finder( 4, 4);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder( 5, 5);
        assertEquals(5, finder.findLeftMostEquivalentPosition());

        finder = finder(6, 6);
        assertEquals(5, finder.findLeftMostEquivalentPosition());

        finder = finder( 7, 7);
        assertEquals(5, finder.findLeftMostEquivalentPosition());

        finder = finder(8, 8);
        assertEquals(5, finder.findLeftMostEquivalentPosition());

        finder = finder(9, 9);
        assertEquals(5, finder.findLeftMostEquivalentPosition());

        finder = finder( 10, 10);
        assertEquals(5, finder.findLeftMostEquivalentPosition());

        finder = finder( 11, 11);
        assertEquals(11, finder.findLeftMostEquivalentPosition());

        finder = finder( 12, 12);
        assertEquals(11, finder.findLeftMostEquivalentPosition());
    }

    @Test
    public void repeatedSingleBaseDeletion()
    {
        genome = new FixedStringGenome("ACCGGGTTTTAAAAAG");
        LeftMostEquivalentDeletionFinder finder = finder(2, 3);
        assertEquals(2, finder.findLeftMostEquivalentPosition());

        finder = finder(4, 5);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(5, 6);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(7, 8);
        assertEquals(7, finder.findLeftMostEquivalentPosition());

        finder = finder(8, 9);
        assertEquals(7, finder.findLeftMostEquivalentPosition());

        finder = finder(9, 10);
        assertEquals(7, finder.findLeftMostEquivalentPosition());

        finder = finder(11, 12);
        assertEquals(11, finder.findLeftMostEquivalentPosition());

        finder = finder(12, 13);
        assertEquals(11, finder.findLeftMostEquivalentPosition());

        finder = finder(13, 14);
        assertEquals(11, finder.findLeftMostEquivalentPosition());

        finder = finder(14, 15);
        assertEquals(11, finder.findLeftMostEquivalentPosition());
    }

    @Test
    public void repeatLength2Deletion()
    {
        genome = new FixedStringGenome("ACCGCGCGTT");
        LeftMostEquivalentDeletionFinder finder = finder(1, 2);
        assertEquals(1, finder.findLeftMostEquivalentPosition());

        finder = finder(3, 4);
        assertEquals(3, finder.findLeftMostEquivalentPosition());

        finder = finder(5, 6);
        assertEquals(3, finder.findLeftMostEquivalentPosition());

        finder = finder(7, 8);
        assertEquals(3, finder.findLeftMostEquivalentPosition());

        finder = finder(9, 10);
        assertEquals(9, finder.findLeftMostEquivalentPosition());
    }

    @Test
    public void doNotSkip()
    {
        genome = new FixedStringGenome("ATGGCGGAAA");
        LeftMostEquivalentDeletionFinder finder = finder(7, 7);
        assertEquals(6, finder.findLeftMostEquivalentPosition());
    }


    @Test
    public void repeatLength3Deletion()
    {
        genome = new FixedStringGenome("AAAGTTGTTGTTGTT");
        LeftMostEquivalentDeletionFinder finder = finder(1, 3);
        assertEquals(1, finder.findLeftMostEquivalentPosition());

        finder = finder(2, 4);
        assertEquals(2, finder.findLeftMostEquivalentPosition());

        finder = finder(3, 5);
        assertEquals(3, finder.findLeftMostEquivalentPosition());

        finder = finder(4, 6);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(5, 7);
        assertEquals(5, finder.findLeftMostEquivalentPosition());

        finder = finder(6, 8);
        assertEquals(6, finder.findLeftMostEquivalentPosition());

        finder = finder(7, 9);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(8, 10);
        assertEquals(5, finder.findLeftMostEquivalentPosition());

        finder = finder(9, 11);
        assertEquals(6, finder.findLeftMostEquivalentPosition());

        finder = finder(10, 12);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(11, 13);
        assertEquals(5, finder.findLeftMostEquivalentPosition());

        finder = finder(12, 14);
        assertEquals(6, finder.findLeftMostEquivalentPosition());

        finder = finder(7, 12);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(8, 13);
        assertEquals(5, finder.findLeftMostEquivalentPosition());
    }

    private LeftMostEquivalentDeletionFinder finder(int start, int end)
    {
        return new LeftMostEquivalentDeletionFinder(genome, "1", start, end);
    }
}

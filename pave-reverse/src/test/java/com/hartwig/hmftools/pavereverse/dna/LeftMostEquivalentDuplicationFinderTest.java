package com.hartwig.hmftools.pavereverse.dna;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.pavereverse.FixedStringGenome;

import org.junit.Test;

public class LeftMostEquivalentDuplicationFinderTest
{
    private RefGenomeInterface genome;
    
    @Test
    public void lookRight()
    {
        // GTCCTCCTCTCCT
        // GTCCTCCT CTC CT
        // >GTCCTCCT CTC CTC CT

        // GTC CTC CTCTCCT
        // <GTC CTC CTC CTCTCCT

        // >GTCCTCCTCTCCTCCT
        // <GTCCTCCTCCTCTCCT

        genome = new FixedStringGenome("GTCCTCCTCTCCT");
        LeftMostEquivalentDuplicationFinder finder = finder( 9, 11);
        assertEquals(9, finder.findLeftMostEquivalentPosition());
    }

    @Test
    public void singleBaseDuplication()
    {
        genome = new FixedStringGenome("CCTGAAAAAAGG");
        LeftMostEquivalentDuplicationFinder finder = finder( 1, 1);
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
    public void repeatedSingleBaseDuplication()
    {
        genome = new FixedStringGenome("ACCGGGTTTTAAAAAAG");
        LeftMostEquivalentDuplicationFinder finder = finder(2, 3);
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

        finder = finder(8, 10);
        assertEquals(7, finder.findLeftMostEquivalentPosition());

        finder = finder(14, 16);
        assertEquals(11, finder.findLeftMostEquivalentPosition());

        finder = finder(13, 16);
        assertEquals(11, finder.findLeftMostEquivalentPosition());

        finder = finder(12, 16);
        assertEquals(11, finder.findLeftMostEquivalentPosition());
    }

    @Test
    public void repeatLength2Duplication()
    {
        genome = new FixedStringGenome("ACCGCGCGTT");
        LeftMostEquivalentDuplicationFinder finder = finder(1, 2);
        assertEquals(1, finder.findLeftMostEquivalentPosition());

        finder = finder(3, 4);
        assertEquals(3, finder.findLeftMostEquivalentPosition());

        finder = finder(5, 6);
        assertEquals(3, finder.findLeftMostEquivalentPosition());

        finder = finder(7, 8);
        assertEquals(3, finder.findLeftMostEquivalentPosition());

        finder = finder(5, 8);
        assertEquals(3, finder.findLeftMostEquivalentPosition());
    }

    @Test
    public void doNotSkip()
    {
        genome = new FixedStringGenome("ATGGCGGAAA");
        LeftMostEquivalentDuplicationFinder finder = finder(7, 7);
        assertEquals(6, finder.findLeftMostEquivalentPosition());
    }

    @Test
    public void repeatLength3Duplication()
    {
        // Note xyz xyz xyz = x yzx yzx yz = xy zxy zxy z
        genome = new FixedStringGenome("AAAGTTGTTGTTGTT");
        LeftMostEquivalentDuplicationFinder finder = finder(1, 3);
        assertEquals(1, finder.findLeftMostEquivalentPosition());

        finder = finder(2, 4);
        assertEquals(2, finder.findLeftMostEquivalentPosition());

        finder = finder(3, 5);
        assertEquals(3, finder.findLeftMostEquivalentPosition());

        finder = finder(4, 6);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(5, 7);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(6, 8);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(7, 9);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(8, 10);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(9, 11);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(10, 12);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(11, 13);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(12, 14);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(7, 12);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(8, 13);
        assertEquals(4, finder.findLeftMostEquivalentPosition());

        finder = finder(7, 15);
        assertEquals(4, finder.findLeftMostEquivalentPosition());
    }

    private LeftMostEquivalentDuplicationFinder finder(int start, int end)
    {
        return new LeftMostEquivalentDuplicationFinder(genome, "1", start, end);
    }
}

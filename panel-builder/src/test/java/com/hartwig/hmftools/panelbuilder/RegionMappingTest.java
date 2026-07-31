package com.hartwig.hmftools.panelbuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThrows;

import java.util.List;
import java.util.OptionalInt;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class RegionMappingTest
{
    // Regions: chr1[100,109] (10b, space [0,10)), chr1[200,214] (15b, space [10,25)), chr2[50,54] (5b, space [25,30)).
    private static final RegionMapping MAPPING = new RegionMapping(List.of(
            new ChrBaseRegion("1", 100, 109),
            new ChrBaseRegion("1", 200, 214),
            new ChrBaseRegion("2", 50, 54)));

    @Test
    public void testLength()
    {
        assertEquals(30, MAPPING.length());
    }

    @Test
    public void testToProbeSpacePosition()
    {
        assertEquals(OptionalInt.of(0), MAPPING.toProbeSpacePosition("1", 100));
        assertEquals(OptionalInt.of(9), MAPPING.toProbeSpacePosition("1", 109));
        assertEquals(OptionalInt.of(10), MAPPING.toProbeSpacePosition("1", 200));
        assertEquals(OptionalInt.of(24), MAPPING.toProbeSpacePosition("1", 214));
        assertEquals(OptionalInt.of(25), MAPPING.toProbeSpacePosition("2", 50));
        // Intron between mapped regions.
        assertFalse(MAPPING.toProbeSpacePosition("1", 150).isPresent());
        // Unmapped chromosome.
        assertFalse(MAPPING.toProbeSpacePosition("3", 100).isPresent());
    }

    @Test
    public void testToGenomePosition()
    {
        assertGenomePosition("1", 100, MAPPING.toGenomePosition(0));
        assertGenomePosition("1", 200, MAPPING.toGenomePosition(10));
        assertGenomePosition("1", 214, MAPPING.toGenomePosition(24));
        assertGenomePosition("2", 54, MAPPING.toGenomePosition(29));
        assertThrows(IllegalArgumentException.class, () -> MAPPING.toGenomePosition(-1));
        assertThrows(IllegalArgumentException.class, () -> MAPPING.toGenomePosition(30));
    }

    @Test
    public void testPositionRoundTrip()
    {
        for(int space = 0; space < MAPPING.length(); ++space)
        {
            BasePosition genome = MAPPING.toGenomePosition(space);
            assertEquals(OptionalInt.of(space), MAPPING.toProbeSpacePosition(genome.Chromosome, genome.Position));
        }
    }

    @Test
    public void testToGenomeRegionsWithinOneRegion()
    {
        assertEquals(List.of(new ChrBaseRegion("1", 100, 109)), MAPPING.toGenomeRegions(0, 10));
        assertEquals(List.of(new ChrBaseRegion("2", 50, 54)), MAPPING.toGenomeRegions(25, 30));
    }

    @Test
    public void testToGenomeRegionsAcrossJunction()
    {
        assertEquals(
                List.of(new ChrBaseRegion("1", 105, 109), new ChrBaseRegion("1", 200, 201)),
                MAPPING.toGenomeRegions(5, 12));
    }

    @Test
    public void testToGenomeRegionsFullSpan()
    {
        assertEquals(
                List.of(
                        new ChrBaseRegion("1", 100, 109),
                        new ChrBaseRegion("1", 200, 214),
                        new ChrBaseRegion("2", 50, 54)),
                MAPPING.toGenomeRegions(0, 30));
    }

    @Test
    public void testToGenomeRegionsInvalidRange()
    {
        assertThrows(IllegalArgumentException.class, () -> MAPPING.toGenomeRegions(-1, 5));
        assertThrows(IllegalArgumentException.class, () -> MAPPING.toGenomeRegions(5, 5));
        assertThrows(IllegalArgumentException.class, () -> MAPPING.toGenomeRegions(0, 31));
    }

    private static void assertGenomePosition(final String chromosome, int position, final BasePosition actual)
    {
        assertEquals(chromosome, actual.Chromosome);
        assertEquals(position, actual.Position);
    }

    @Test
    public void testConstructorInvalid()
    {
        assertThrows(IllegalArgumentException.class, () -> new RegionMapping(List.of()));
        // Overlapping regions.
        assertThrows(
                IllegalArgumentException.class, () -> new RegionMapping(List.of(
                        new ChrBaseRegion("1", 100, 109),
                        new ChrBaseRegion("1", 105, 114))));
        // Not in ascending order.
        assertThrows(
                IllegalArgumentException.class, () -> new RegionMapping(List.of(
                        new ChrBaseRegion("1", 200, 214),
                        new ChrBaseRegion("1", 100, 109))));
    }
}

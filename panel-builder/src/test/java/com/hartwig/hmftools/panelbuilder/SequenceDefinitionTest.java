package com.hartwig.hmftools.panelbuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class SequenceDefinitionTest
{
    public static final ChrBaseRegion REGION1 = new ChrBaseRegion("1", 1, 10);
    public static final ChrBaseRegion REGION2 = new ChrBaseRegion("2", 21, 30);
    public static final String INSERT = "ACGT";

    @Test
    public void testFactoriesValid()
    {
        SequenceDefinition.singleRegion(REGION1);
        SequenceDefinition.forwardSgl(REGION1, INSERT);
        SequenceDefinition.reverseSgl(INSERT, REGION2);
        SequenceDefinition.variant(REGION1, Orientation.REVERSE, "", REGION2, Orientation.FORWARD);
        SequenceDefinition.variant(REGION1, Orientation.REVERSE, INSERT, REGION2, Orientation.FORWARD);
    }

    @Test
    public void testConstructorInvalid()
    {
        // No segments.
        assertThrows(IllegalArgumentException.class, () -> new SequenceDefinition(List.of()));

        // Insert only (no region).
        assertThrows(IllegalArgumentException.class, () -> new SequenceDefinition(List.of(new InsertSeqSegment(INSERT))));

        // Consecutive inserts should be a single insert segment.
        assertThrows(
                IllegalArgumentException.class, () -> new SequenceDefinition(List.of(
                        new RefSegment(REGION1, Orientation.FORWARD),
                        new InsertSeqSegment("AC"),
                        new InsertSeqSegment("GT"))));

        // Adjacent regions with the same orientation should be a single region.
        assertThrows(
                IllegalArgumentException.class, () -> new SequenceDefinition(List.of(
                        new RefSegment(new ChrBaseRegion("1", 1, 10), Orientation.FORWARD),
                        new RefSegment(new ChrBaseRegion("1", 11, 20), Orientation.FORWARD))));
        assertThrows(
                IllegalArgumentException.class, () -> new SequenceDefinition(List.of(
                        new RefSegment(new ChrBaseRegion("1", 11, 20), Orientation.REVERSE),
                        new RefSegment(new ChrBaseRegion("1", 1, 10), Orientation.REVERSE))));
    }

    @Test
    public void testAdjacentRegionsAllowed()
    {
        // Adjacent in the genome but different orientation, so a genuine junction.
        new SequenceDefinition(List.of(
                new RefSegment(new ChrBaseRegion("1", 1, 10), Orientation.FORWARD),
                new RefSegment(new ChrBaseRegion("1", 11, 20), Orientation.REVERSE)));
        // Adjacent positions but different chromosomes.
        new SequenceDefinition(List.of(
                new RefSegment(new ChrBaseRegion("1", 1, 10), Orientation.FORWARD),
                new RefSegment(new ChrBaseRegion("2", 11, 20), Orientation.FORWARD)));
    }

    @Test
    public void testSegmentValidation()
    {
        assertThrows(IllegalArgumentException.class, () -> new InsertSeqSegment(""));
        assertThrows(IllegalArgumentException.class, () -> new InsertSeqSegment("ACGTN"));
        assertThrows(IllegalArgumentException.class, () -> new RefSegment(new ChrBaseRegion("1", 10, 1), Orientation.FORWARD));
    }

    @Test
    public void testSingleRegion()
    {
        SequenceDefinition actual = SequenceDefinition.singleRegion(REGION1);
        assertTrue(actual.isSingleRegion());
        assertFalse(actual.isMultiRegion());
        assertEquals(REGION1, actual.singleRegion());
        assertEquals(REGION1, actual.singleRegionOrNull());
        assertEquals(List.of(REGION1), actual.regions());
        assertEquals(10, actual.baseLength());
    }

    @Test
    public void testStructuralVariant()
    {
        SequenceDefinition actual = SequenceDefinition.variant(REGION1, Orientation.REVERSE, INSERT, REGION2, Orientation.FORWARD);
        assertFalse(actual.isSingleRegion());
        assertTrue(actual.isMultiRegion());
        assertThrows(IllegalArgumentException.class, actual::singleRegion);
        assertNull(actual.singleRegionOrNull());
        assertEquals(List.of(REGION1, REGION2), actual.regions());
        assertEquals(24, actual.baseLength());
    }

    @Test
    public void testSplicedMultiRegion()
    {
        ChrBaseRegion region3 = new ChrBaseRegion("1", 41, 45);
        SequenceDefinition actual = new SequenceDefinition(List.of(
                new RefSegment(REGION1, Orientation.FORWARD),
                new RefSegment(REGION2, Orientation.FORWARD),
                new RefSegment(region3, Orientation.FORWARD)));
        assertFalse(actual.isSingleRegion());
        assertTrue(actual.isMultiRegion());
        assertThrows(IllegalArgumentException.class, actual::singleRegion);
        assertNull(actual.singleRegionOrNull());
        assertEquals(List.of(REGION1, REGION2, region3), actual.regions());
        assertEquals(25, actual.baseLength());
    }

    @Test
    public void testSplicedFactory()
    {
        ChrBaseRegion region3 = new ChrBaseRegion("1", 41, 45);
        SequenceDefinition actual = SequenceDefinition.spliced(List.of(REGION1, REGION2, region3));
        assertTrue(actual.isMultiRegion());
        assertEquals(List.of(REGION1, REGION2, region3), actual.regions());
        assertEquals(25, actual.baseLength());
        // Single region collapses to a plain single-region definition.
        assertTrue(SequenceDefinition.spliced(List.of(REGION1)).isSingleRegion());
    }

    @Test
    public void testCompareTo()
    {
        SequenceDefinition a = SequenceDefinition.singleRegion(new ChrBaseRegion("1", 1, 10));
        SequenceDefinition b = SequenceDefinition.singleRegion(new ChrBaseRegion("1", 20, 30));
        // Ordered by region position.
        assertTrue(a.compareTo(b) < 0);
        assertTrue(b.compareTo(a) > 0);
        assertEquals(0, a.compareTo(SequenceDefinition.singleRegion(new ChrBaseRegion("1", 1, 10))));

        // Shorter tuple orders first when it is a prefix of the other.
        SequenceDefinition sgl = SequenceDefinition.forwardSgl(new ChrBaseRegion("1", 1, 10), "AC");
        assertTrue(a.compareTo(sgl) < 0);

        // Cross-type second segment: a region (rank 0) orders before an insert (rank 1).
        SequenceDefinition twoRegions = SequenceDefinition.variant(
                new ChrBaseRegion("1", 1, 10), Orientation.FORWARD, "", new ChrBaseRegion("2", 1, 10), Orientation.FORWARD);
        assertTrue(twoRegions.compareTo(sgl) < 0);

        // Same region, orientation breaks the tie (FORWARD before REVERSE).
        SequenceDefinition reverse = new SequenceDefinition(List.of(new RefSegment(new ChrBaseRegion("1", 1, 10), Orientation.REVERSE)));
        assertTrue(a.compareTo(reverse) < 0);
    }

    @Test
    public void testForwardSgl()
    {
        SequenceDefinition actual = SequenceDefinition.forwardSgl(REGION1, INSERT);
        assertFalse(actual.isSingleRegion());
        assertFalse(actual.isMultiRegion());
        assertThrows(IllegalArgumentException.class, actual::singleRegion);
        assertNull(actual.singleRegionOrNull());
        assertEquals(List.of(REGION1), actual.regions());
        assertEquals(14, actual.baseLength());
    }

    @Test
    public void testReverseSgl()
    {
        SequenceDefinition actual = SequenceDefinition.reverseSgl(INSERT, REGION2);
        assertFalse(actual.isSingleRegion());
        assertFalse(actual.isMultiRegion());
        assertThrows(IllegalArgumentException.class, actual::singleRegion);
        assertNull(actual.singleRegionOrNull());
        assertEquals(List.of(REGION2), actual.regions());
        assertEquals(14, actual.baseLength());
    }
}

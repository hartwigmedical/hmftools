package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import org.junit.Test;

public class VariantTest
{
    @Test
    public void testPointMutation()
    {
        // SNVs and MNVs
        VariantData var = new VariantData(CHR_1, 100, "A", "G");

        assertTrue(var.isBaseChange());
        assertEquals(100, var.EndPosition);
        assertEquals(1, var.altPositions().size());
        assertTrue(var.altBasesAbove(99));
        assertFalse(var.altBasesAbove(100));
        assertTrue(var.altBasesBelow(101));
        assertFalse(var.altBasesBelow(100));

        assertTrue(var.altPositionsOverlap(90, 100));
        assertFalse(var.altPositionsOverlap(90, 99));
        assertTrue(var.altPositionsOverlap(100, 110));
        assertFalse(var.altPositionsOverlap(101, 110));

        assertTrue(var.altPositionsWithin(100, 110));
        assertTrue(var.altPositionsWithin(90, 100));

        var = new VariantData(CHR_1, 100, "AGA", "GCG");

        assertTrue(var.isBaseChange());
        assertEquals(102, var.EndPosition);
        assertEquals(3, var.altPositions().size());
        assertTrue(var.altBasesAbove(99));
        assertFalse(var.altBasesAbove(100));
        assertFalse(var.altBasesAbove(102));
        assertTrue(var.altBasesBelow(103));
        assertFalse(var.altBasesBelow(102));
        assertFalse(var.altBasesBelow(100));

        assertTrue(var.altPositionsOverlap(90, 100));
        assertFalse(var.altPositionsOverlap(103, 110));

        assertTrue(var.altPositionsWithin(100, 110));
        assertFalse(var.altPositionsWithin(102, 110));
        assertTrue(var.altPositionsWithin(90, 102));
        assertFalse(var.altPositionsWithin(90, 101));
    }

    @Test
    public void testDeletion()
    {
        VariantData var = new VariantData(CHR_1, 100, "ACGT", "A");

        assertTrue(var.isDeletion());
        assertEquals(104, var.EndPosition);
        assertEquals(3, var.altPositions().size());
        assertTrue(var.altBasesAbove(99));
        assertTrue(var.altBasesAbove(100));
        assertFalse(var.altBasesAbove(101));
        assertTrue(var.altBasesBelow(104));
        assertFalse(var.altBasesBelow(103));

        assertFalse(var.altPositionsOverlap(90, 100));
        assertTrue(var.altPositionsOverlap(101, 110));
        assertFalse(var.altPositionsOverlap(104, 110));
        assertTrue(var.altPositionsOverlap(103, 110));

        assertTrue(var.altPositionsWithin(100, 110));
        assertTrue(var.altPositionsWithin(101, 110));
        assertFalse(var.altPositionsWithin(102, 110));
        assertTrue(var.altPositionsWithin(90, 104));
        assertTrue(var.altPositionsWithin(90, 103));
        assertFalse(var.altPositionsWithin(90, 102));
    }

    @Test
    public void testInsertion()
    {
        VariantData var = new VariantData(CHR_1, 100, "A", "AAAAA"); // 4 bases

        assertTrue(var.isInsert());
        assertEquals(101, var.EndPosition);
        assertTrue(var.altPositions().isEmpty());

        assertTrue(var.altBasesAbove(100));
        assertFalse(var.altBasesAbove(101));
        assertTrue(var.altBasesBelow(101));
        assertFalse(var.altBasesBelow(100));

        assertFalse(var.altPositionsOverlap(90, 100));
        assertTrue(var.altPositionsOverlap(90, 101));
        assertFalse(var.altPositionsOverlap(101, 110));
        assertTrue(var.altPositionsOverlap(100, 110));

        assertTrue(var.altPositionsWithin(100, 110));
        assertFalse(var.altPositionsWithin(101, 110));
        assertTrue(var.altPositionsWithin(90, 101));
        assertFalse(var.altPositionsWithin(90, 100));
    }
}

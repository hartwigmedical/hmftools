package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.selfAlignment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

// The placement invariant: a record has candidates and a chosen one among them, or neither.
public class LiftedRecordTest
{
    private static LiftedRecord resolved(final int primaryIndex, final List<LiftedAlignment> alignments)
    {
        return new LiftedRecord(60, 1, "", primaryIndex, alignments);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testResolvedWithoutPlacementRejected()
    {
        resolved(LiftedRecord.NO_PRIMARY, List.of(selfAlignment(CHR_1, 1000, "50M")));
    }

    @Test(expected = IllegalArgumentException.class)
    public void testResolvedWithoutAlignmentsRejected()
    {
        resolved(0, List.of());
    }

    @Test(expected = IllegalArgumentException.class)
    public void testPrimaryIndexPastEndRejected()
    {
        resolved(2, List.of(selfAlignment(CHR_1, 1000, "50M")));
    }

    @Test(expected = IllegalArgumentException.class)
    public void testAlignmentsWithoutAChosenPlacementRejected()
    {
        new LiftedRecord(0, 0, "", LiftedRecord.NO_PRIMARY, List.of(selfAlignment(CHR_1, 1000, "50M")));
    }

    @Test
    public void testPlacementReadFromTheIndexedAlignment()
    {
        LiftedRecord liftedRecord = resolved(
                1, List.of(selfAlignment(CHR_1, 1000, "50M"), selfAlignment(CHR_1, 5000, "20M100N30M")));

        assertTrue(liftedRecord.hasPlacement());
        assertTrue(liftedRecord.swapped());
        assertEquals(5000, liftedRecord.finalPos());
        assertEquals("20M100N30M", liftedRecord.finalCigar());
        assertTrue(liftedRecord.hasNCigar());
    }

    @Test
    public void testUnmappedHasNoPlacement()
    {
        LiftedRecord liftedRecord = LiftedRecord.unmapped("over_cap_unmapped");

        assertFalse(liftedRecord.hasPlacement());
        assertFalse(liftedRecord.swapped());
        assertEquals("over_cap_unmapped", liftedRecord.notes());
    }
}

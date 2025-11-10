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
    public void testConstructorValid()
    {
        // Single region.
        new SequenceDefinition(REGION1, null, "", null, null);
        // SGL.
        new SequenceDefinition(REGION1, Orientation.FORWARD, INSERT, null, null);
        new SequenceDefinition(null, null, INSERT, REGION2, Orientation.FORWARD);
        // Generic variant.
        new SequenceDefinition(REGION1, Orientation.REVERSE, "", REGION2, Orientation.FORWARD);
        new SequenceDefinition(REGION1, Orientation.REVERSE, INSERT, REGION2, Orientation.FORWARD);
    }

    @Test
    public void testConstructorInvalid()
    {
        // Single region but wrong format.
        assertThrows(IllegalArgumentException.class, () -> new SequenceDefinition(REGION1, Orientation.FORWARD, "", null, null));
        assertThrows(IllegalArgumentException.class, () -> new SequenceDefinition(null, null, "", REGION2, null));
        assertThrows(IllegalArgumentException.class, () -> new SequenceDefinition(null, null, "", REGION2, Orientation.FORWARD));

        // Orientation missing.
        assertThrows(IllegalArgumentException.class, () -> new SequenceDefinition(REGION1, null, INSERT, null, null));
        assertThrows(IllegalArgumentException.class, () -> new SequenceDefinition(null, null, INSERT, REGION2, null));

        // Adjacent regions with no insert.
        assertThrows(IllegalArgumentException.class, () -> new SequenceDefinition(new ChrBaseRegion("1", 1, 10), Orientation.FORWARD, "", new ChrBaseRegion("1", 11, 20), Orientation.FORWARD));

        // Insert only.
        assertThrows(IllegalArgumentException.class, () -> new SequenceDefinition(null, null, INSERT, null, null));

        // Everything empty.
        assertThrows(IllegalArgumentException.class, () -> new SequenceDefinition(null, null, "", null, null));
    }

    @Test
    public void testSingleRegion()
    {
        SequenceDefinition actual = SequenceDefinition.singleRegion(REGION1);
        SequenceDefinition expected = new SequenceDefinition(REGION1, null, "", null, null);
        assertEquals(expected, actual);
        assertTrue(actual.isSingleRegion());
        assertEquals(REGION1, actual.singleRegion());
        assertEquals(REGION1, actual.singleRegionOrNull());
        assertEquals(List.of(REGION1), actual.regions());
        assertEquals(10, actual.baseLength());
    }

    @Test
    public void testStructuralVariant()
    {
        SequenceDefinition actual = new SequenceDefinition(REGION1, Orientation.REVERSE, INSERT, REGION2, Orientation.FORWARD);
        assertFalse(actual.isSingleRegion());
        assertThrows(IllegalArgumentException.class, actual::singleRegion);
        assertNull(actual.singleRegionOrNull());
        assertEquals(List.of(REGION1, REGION2), actual.regions());
        assertEquals(24, actual.baseLength());
    }

    @Test
    public void testForwardSgl()
    {
        SequenceDefinition actual =
                SequenceDefinition.forwardSgl(REGION1, INSERT);
        SequenceDefinition expected = new SequenceDefinition(REGION1, Orientation.FORWARD, INSERT, null, null);
        assertEquals(expected, actual);
        assertFalse(actual.isSingleRegion());
        assertThrows(IllegalArgumentException.class, actual::singleRegion);
        assertNull(actual.singleRegionOrNull());
        assertEquals(List.of(REGION1), actual.regions());
        assertEquals(14, actual.baseLength());
    }

    @Test
    public void testReverseSgl()
    {
        SequenceDefinition actual =
                SequenceDefinition.reverseSgl(INSERT, REGION2);
        SequenceDefinition expected = new SequenceDefinition(null, null, INSERT, REGION2, Orientation.FORWARD);
        assertEquals(expected, actual);
        assertFalse(actual.isSingleRegion());
        assertThrows(IllegalArgumentException.class, actual::singleRegion);
        assertNull(actual.singleRegionOrNull());
        assertEquals(List.of(REGION2), actual.regions());
        assertEquals(14, actual.baseLength());
    }
}

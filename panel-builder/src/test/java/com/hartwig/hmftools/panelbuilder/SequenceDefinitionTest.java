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
    public void testConstructor()
    {
        new SequenceDefinition(REGION1, null, null, null, null);
        new SequenceDefinition(REGION1, null, INSERT, null, null);
        new SequenceDefinition(null, null, INSERT, REGION2, null);
        new SequenceDefinition(REGION1, null, INSERT, REGION2, null);

        assertThrows(IllegalArgumentException.class, () -> new SequenceDefinition(null, null, null, null, null));
        assertThrows(IllegalArgumentException.class, () -> new SequenceDefinition(REGION1, null, null, REGION2, null));
        assertThrows(IllegalArgumentException.class, () -> new SequenceDefinition(null, null, null, REGION2, null));
        assertThrows(IllegalArgumentException.class, () -> new SequenceDefinition(null, null, INSERT, null, null));
    }

    @Test
    public void testExactRegion()
    {
        SequenceDefinition actual = SequenceDefinition.exactRegion(REGION1);
        SequenceDefinition expected = new SequenceDefinition(REGION1, null, null, null, null);
        assertEquals(expected, actual);
        assertTrue(actual.isExactRegion());
        assertEquals(REGION1, actual.exactRegion());
        assertEquals(REGION1, actual.exactRegionOrNull());
        assertEquals(List.of(REGION1), actual.regions());
        assertEquals(10, actual.baseLength());
    }

    @Test
    public void testSimpleMutation()
    {
        SequenceDefinition actual = SequenceDefinition.simpleMutation(REGION1, INSERT, REGION2);
        SequenceDefinition expected = new SequenceDefinition(REGION1, Orientation.FORWARD, INSERT, REGION2, Orientation.FORWARD);
        assertEquals(expected, actual);
        assertFalse(actual.isExactRegion());
        assertThrows(IllegalArgumentException.class, actual::exactRegion);
        assertNull(actual.exactRegionOrNull());
        assertEquals(List.of(REGION1, REGION2), actual.regions());
        assertEquals(24, actual.baseLength());
    }

    @Test
    public void testStructuralVariant()
    {
        SequenceDefinition actual =
                SequenceDefinition.structuralVariant(REGION1, Orientation.REVERSE, INSERT, REGION2, Orientation.FORWARD);
        SequenceDefinition expected = new SequenceDefinition(REGION1, Orientation.REVERSE, INSERT, REGION2, Orientation.FORWARD);
        assertEquals(expected, actual);
        assertFalse(actual.isExactRegion());
        assertThrows(IllegalArgumentException.class, actual::exactRegion);
        assertNull(actual.exactRegionOrNull());
        assertEquals(List.of(REGION1, REGION2), actual.regions());
        assertEquals(24, actual.baseLength());
    }
}

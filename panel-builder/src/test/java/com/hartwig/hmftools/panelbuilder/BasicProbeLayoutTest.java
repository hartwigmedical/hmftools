package com.hartwig.hmftools.panelbuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertThrows;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class BasicProbeLayoutTest
{
    private static final ChrBaseRegion REGION1 = new ChrBaseRegion("1", 1, 10);
    private static final ChrBaseRegion REGION2 = new ChrBaseRegion("2", 21, 30);
    private static final String INSERT = "ACGT";

    @Test
    public void testSingleRegion()
    {
        BasicProbeLayout layout = BasicProbeLayout.from(SequenceDefinition.singleRegion(REGION1));
        assertEquals(REGION1, layout.startRegion());
        assertEquals(Orientation.FORWARD, layout.startOrientation());
        assertEquals("", layout.insertSequence());
        assertNull(layout.endRegion());
        assertNull(layout.endOrientation());
    }

    @Test
    public void testForwardSgl()
    {
        BasicProbeLayout layout = BasicProbeLayout.from(SequenceDefinition.forwardSgl(REGION1, INSERT));
        assertEquals(REGION1, layout.startRegion());
        assertEquals(INSERT, layout.insertSequence());
        assertNull(layout.endRegion());
    }

    @Test
    public void testReverseSgl()
    {
        BasicProbeLayout layout = BasicProbeLayout.from(SequenceDefinition.reverseSgl(INSERT, REGION2));
        assertNull(layout.startRegion());
        assertEquals(INSERT, layout.insertSequence());
        assertEquals(REGION2, layout.endRegion());
    }

    @Test
    public void testVariantNoInsert()
    {
        BasicProbeLayout layout = BasicProbeLayout.from(
                SequenceDefinition.variant(REGION1, Orientation.REVERSE, "", REGION2, Orientation.FORWARD));
        assertEquals(REGION1, layout.startRegion());
        assertEquals(Orientation.REVERSE, layout.startOrientation());
        assertEquals("", layout.insertSequence());
        assertEquals(REGION2, layout.endRegion());
        assertEquals(Orientation.FORWARD, layout.endOrientation());
    }

    @Test
    public void testVariantWithInsert()
    {
        BasicProbeLayout layout = BasicProbeLayout.from(
                SequenceDefinition.variant(REGION1, Orientation.FORWARD, INSERT, REGION2, Orientation.FORWARD));
        assertEquals(REGION1, layout.startRegion());
        assertEquals(INSERT, layout.insertSequence());
        assertEquals(REGION2, layout.endRegion());
    }

    @Test
    public void testRejectsMoreThanTwoRegions()
    {
        SequenceDefinition spliced = new SequenceDefinition(List.of(
                new RefSegment(REGION1, Orientation.FORWARD),
                new RefSegment(REGION2, Orientation.FORWARD),
                new RefSegment(new ChrBaseRegion("3", 41, 50), Orientation.FORWARD)));
        assertThrows(IllegalArgumentException.class, () -> BasicProbeLayout.from(spliced));
    }
}

package com.hartwig.hmftools.purple.tools;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;

import org.junit.Test;

public class GeneAmplificationTest
{
    @Test(expected = IllegalArgumentException.class)
    public void testEmptyExonsThrowsException()
    {
        new GeneAmplification(Lists.newArrayList(), 100, 200);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testInvalidRangeThrowsException()
    {
        List<ExonData> exons = Lists.newArrayList(new ExonData(1, 100, 150, 1, 0, 0));
        new GeneAmplification(exons, 200, 100);
    }

    @Test
    public void allExonsAmplified()
    {
        List<ExonData> exons = Lists.newArrayList(
                new ExonData(1, 100, 150, 1, 0, 0),
                new ExonData(1, 200, 250, 2, 0, 0)
        );
        GeneAmplification ga = new GeneAmplification(exons, 50, 300);
        assertTrue(ga.isCompleteAmplification());
        assertEquals(2, ga.numberOfAffectedExons());
        assertTrue(ga.isOfInterest());
        assertFalse(ga.isHeadAmplification());
        assertFalse(ga.isTailAmplification());
    }

    @Test
    public void innerExonAmplified()
    {
        List<ExonData> exons = Lists.newArrayList(
                new ExonData(1, 100, 150, 1, 0, 0),
                new ExonData(1, 200, 250, 2, 0, 0),
                new ExonData(1, 300, 350, 3, 0, 0)
        );
        GeneAmplification ga = new GeneAmplification(exons, 180, 270);
        assertTrue(ga.isOfInterest());
        assertFalse(ga.isCompleteAmplification());
        assertEquals(1, ga.numberOfAffectedExons());
        assertFalse(ga.isHeadAmplification());
        assertFalse(ga.isTailAmplification());
    }

    @Test
    public void noExonsAmplified()
    {
        List<ExonData> exons = Lists.newArrayList(
                new ExonData(1, 100, 150, 1, 0, 0),
                new ExonData(1, 200, 250, 2, 0, 0)
        );
        GeneAmplification ga = new GeneAmplification(exons, 160, 190);
        assertFalse(ga.isOfInterest());
        assertFalse(ga.isCompleteAmplification());
        assertEquals(0, ga.numberOfAffectedExons());
        assertFalse(ga.isHeadAmplification());
        assertFalse(ga.isTailAmplification());
    }

    @Test
    public void firstExonOnlyAmplified()
    {
        List<ExonData> exons = Lists.newArrayList(
                new ExonData(1, 100, 150, 1, 0, 0),
                new ExonData(1, 200, 250, 2, 0, 0),
                new ExonData(1, 300, 350, 3, 0, 0)
        );
        GeneAmplification ga = new GeneAmplification(exons, 50, 170);
        assertFalse(ga.isOfInterest());
        assertFalse(ga.isCompleteAmplification());
        assertEquals(1, ga.numberOfAffectedExons());
        assertTrue(ga.isHeadAmplification());
        assertFalse(ga.isTailAmplification());
    }

    @Test
    public void lastExonOnlyAmplified()
    {
        List<ExonData> exons = Lists.newArrayList(
                new ExonData(1, 100, 150, 1, 0, 0),
                new ExonData(1, 200, 250, 2, 0, 0),
                new ExonData(1, 300, 350, 3, 0, 0)
        );
        GeneAmplification ga = new GeneAmplification(exons, 280, 400);
        assertFalse(ga.isOfInterest());
        assertFalse(ga.isCompleteAmplification());
        assertEquals(1, ga.numberOfAffectedExons());
        assertFalse(ga.isHeadAmplification());
        assertTrue(ga.isTailAmplification());
    }

    @Test
    public void headAmplification()
    {
        List<ExonData> exons = Lists.newArrayList(
                new ExonData(1, 100, 150, 1, 0, 0),
                new ExonData(1, 200, 250, 2, 0, 0),
                new ExonData(1, 300, 350, 3, 0, 0)
        );
        GeneAmplification ga = new GeneAmplification(exons, 50, 270);
        assertFalse(ga.isOfInterest());
        assertFalse(ga.isCompleteAmplification());
        assertEquals(2, ga.numberOfAffectedExons());
        assertTrue(ga.isHeadAmplification());
        assertFalse(ga.isTailAmplification());
    }

    @Test
    public void tailAmplification()
    {
        List<ExonData> exons = Lists.newArrayList(
                new ExonData(1, 100, 150, 1, 0, 0),
                new ExonData(1, 200, 250, 2, 0, 0),
                new ExonData(1, 300, 350, 3, 0, 0)
        );
        GeneAmplification ga = new GeneAmplification(exons, 180, 400);
        assertFalse(ga.isOfInterest());
        assertFalse(ga.isCompleteAmplification());
        assertEquals(2, ga.numberOfAffectedExons());
        assertFalse(ga.isHeadAmplification());
        assertTrue(ga.isTailAmplification());
    }
}

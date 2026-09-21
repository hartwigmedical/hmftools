package com.hartwig.hmftools.virusdetect.common;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.util.List;

import com.hartwig.hmftools.virusdetect.app.SummaryStats;

import org.junit.Test;

public class SummaryStatsTest
{
    private static final double EPSILON = 1e-9;

    @Test
    public void testMeanMinMaxAndPercentiles()
    {
        SummaryStats stats = SummaryStats.from(new int[] { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 });

        assertEquals(55.0, stats.mean(), EPSILON);
        assertEquals(10.0, stats.min(), EPSILON);
        assertEquals(100.0, stats.max(), EPSILON);
        assertEquals(10.0, stats.p5(), EPSILON);    // ceil(0.05*10)-1 = 0
        assertEquals(30.0, stats.p25(), EPSILON);   // ceil(0.25*10)-1 = 2
        assertEquals(50.0, stats.p50(), EPSILON);   // ceil(0.50*10)-1 = 4
        assertEquals(80.0, stats.p75(), EPSILON);   // ceil(0.75*10)-1 = 7
        assertEquals(100.0, stats.p95(), EPSILON);  // ceil(0.95*10)-1 = 9
    }

    @Test
    public void testSingleValue()
    {
        SummaryStats stats = SummaryStats.from(new int[] { 7 });

        assertEquals(7.0, stats.mean(), EPSILON);
        assertEquals(7.0, stats.min(), EPSILON);
        assertEquals(7.0, stats.max(), EPSILON);
        assertEquals(7.0, stats.p5(), EPSILON);
        assertEquals(7.0, stats.p50(), EPSILON);
        assertEquals(7.0, stats.p95(), EPSILON);
    }

    @Test
    public void testUnsortedInputIsSorted()
    {
        SummaryStats stats = SummaryStats.from(new int[] { 30, 10, 20 });

        assertEquals(20.0, stats.mean(), EPSILON);
        assertEquals(10.0, stats.min(), EPSILON);
        assertEquals(30.0, stats.max(), EPSILON);
        assertEquals(20.0, stats.p50(), EPSILON);
    }

    @Test
    public void testFromListMatchesArray()
    {
        SummaryStats stats = SummaryStats.from(List.of(3, 1, 2));

        assertEquals(2.0, stats.mean(), EPSILON);
        assertEquals(1.0, stats.min(), EPSILON);
        assertEquals(3.0, stats.max(), EPSILON);
        assertEquals(2.0, stats.p50(), EPSILON);
    }

    @Test
    public void testDoesNotMutateInput()
    {
        int[] input = { 30, 10, 20 };
        SummaryStats.from(input);
        assertArrayEquals(new int[] { 30, 10, 20 }, input);
    }

    @Test
    public void testThrowsOnEmpty()
    {
        assertThrows(IllegalArgumentException.class, () -> SummaryStats.from(new int[0]));
    }
}

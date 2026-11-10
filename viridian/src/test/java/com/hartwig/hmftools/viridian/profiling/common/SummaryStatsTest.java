package com.hartwig.hmftools.viridian.profiling.common;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.util.List;

import org.junit.Test;

public class SummaryStatsTest
{
    @Test
    public void testFrom()
    {
        // Nearest-rank percentiles over the sorted values: index ceil(percent * 10) - 1, so p5 at 0, p25 at 2, p50 at 4,
        // p75 at 7 and p95 at 9.

        SummaryStats stats = SummaryStats.from(new int[] { 60, 10, 90, 30, 100, 20, 50, 80, 40, 70 });

        assertEquals(new SummaryStats(55.0, 10.0, 10.0, 30.0, 50.0, 80.0, 100.0, 100.0), stats);
    }

    @Test
    public void testFromList()
    {
        assertEquals(SummaryStats.from(new int[] { 30, 10, 20 }), SummaryStats.from(List.of(30, 10, 20)));
    }

    @Test
    public void testFromInputMutation()
    {
        int[] input = { 30, 10, 20 };
        SummaryStats.from(input);
        assertArrayEquals(new int[] { 30, 10, 20 }, input);
    }

    @Test
    public void testFromEmpty()
    {
        assertThrows(IllegalArgumentException.class, () -> SummaryStats.from(new int[0]));
    }
}

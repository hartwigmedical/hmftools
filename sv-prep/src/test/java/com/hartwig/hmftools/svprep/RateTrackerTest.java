package com.hartwig.hmftools.svprep;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import org.junit.Test;

public class RateTrackerTest
{
    @Test
    public void testRateTracker()
    {
        ReadRateTracker readRateTracker = new ReadRateTracker(10, 1, 5);

        assertTrue(readRateTracker.handleRead(1));
        assertTrue(readRateTracker.handleRead(2));
        assertTrue(readRateTracker.handleRead(3));

        // reset at next segment
        assertTrue(readRateTracker.handleRead(11));
        assertEquals(1, readRateTracker.readCount());
        assertFalse(readRateTracker.isRateLimited());

        // exceed the limit in this segment
        for(int i = 0; i < 14; ++i)
        {
            assertTrue(readRateTracker.handleRead(12));
        }

        assertFalse(readRateTracker.isRateLimited());

        assertFalse(readRateTracker.handleRead(21));
        assertTrue(readRateTracker.isRateLimited());

        assertFalse(readRateTracker.handleRead(21));
        assertTrue(readRateTracker.handleRead(21));

        assertFalse(readRateTracker.handleRead(24));
        assertFalse(readRateTracker.handleRead(25));
        assertTrue(readRateTracker.handleRead(26));

        // resets but still above
        assertFalse(readRateTracker.handleRead(31));
        assertTrue(readRateTracker.isRateLimited());

        // resets and below
        assertTrue(readRateTracker.handleRead(42));
        assertFalse(readRateTracker.isRateLimited());
    }
}

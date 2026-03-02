package com.hartwig.hmftools.amber.contamination;

import static java.util.List.of;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Assert;
import org.junit.Test;

public class VafConsistencyCheckResultTest
{
    @Test
    public void proportionTest()
    {
        assertEquals(0.001, new VafConsistencyCheckResult(0.1, 10, 10000, of()).proportion(), 0.0001);
        assertEquals(0.01, new VafConsistencyCheckResult(0.1, 1, 100, of()).proportion(), 0.0001);
        assertEquals(0.1, new VafConsistencyCheckResult(0.1, 10, 100, of()).proportion(), 0.0001);
    }
}

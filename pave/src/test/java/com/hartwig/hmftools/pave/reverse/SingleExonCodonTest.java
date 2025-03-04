package com.hartwig.hmftools.pave.reverse;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class SingleExonCodonTest extends ReversePaveTestBase
{
    @Test
    public void strandLocationOfChangeTest()
    {
        assertEquals(1000, sec(1000, "AGA", true).forwardStrandLocationOfChange(0));
        assertEquals(1001, sec(1000, "AGA", true).forwardStrandLocationOfChange(1));
        assertEquals(1002, sec(1000, "AGA", true).forwardStrandLocationOfChange(2));
    }
}

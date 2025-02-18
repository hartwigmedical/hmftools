package com.hartwig.hmftools.pave.transval;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class SingleExonCodonTest extends TransvalTestBase
{
    @Test
    public void strandLocationOfChangeTest()
    {
        assertEquals(1000, sec(1000, "AGA").strandLocationOfChange(0));
        assertEquals(1001, sec(1000, "AGA").strandLocationOfChange(1));
        assertEquals(1002, sec(1000, "AGA").strandLocationOfChange(2));
    }
}

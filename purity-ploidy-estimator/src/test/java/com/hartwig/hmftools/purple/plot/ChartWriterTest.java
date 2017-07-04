package com.hartwig.hmftools.purple.plot;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ChartWriterTest {

    @Test
    public void testSubtitle() {
        assertEquals("COLO829 P:100% NF:0.64", ChartWriter.subtitle("COLO829", 1.0, 0.64));
        assertEquals("CPCT012345 P:83% NF:1.10", ChartWriter.subtitle("CPCT012345", 0.83, 1.1));
        assertEquals("CPCT012345 P:83% NF:1.10", ChartWriter.subtitle("CPCT012345", 0.83, 1.101));
    }

}

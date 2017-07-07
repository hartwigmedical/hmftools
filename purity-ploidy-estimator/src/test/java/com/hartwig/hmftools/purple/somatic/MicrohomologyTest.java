package com.hartwig.hmftools.purple.somatic;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class MicrohomologyTest {

    @Test
    public void testDelete() {
        assertEquals("TATC", Microhomology.microhomology("TATCCACTTGAT", "T", "TATCCACTTGATATCATTAT"));
    }
}

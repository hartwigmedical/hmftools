package com.hartwig.hmftools.common.codon;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class CodonsTest {

    @Test
    public void testAsCodon() {
        assertEquals('Y', Codons.asCondon("TAT"));
    }

    @Test
    public void testAsCodonString() {
        assertEquals("Y", Codons.asCodonString("TAT"));
        assertEquals("Y", Codons.asCodonString("TATG"));
        assertEquals("Y", Codons.asCodonString("TATGA"));
        assertEquals("YD", Codons.asCodonString("TATGAT"));
    }

}

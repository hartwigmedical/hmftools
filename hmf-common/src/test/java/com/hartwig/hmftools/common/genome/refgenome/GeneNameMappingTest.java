package com.hartwig.hmftools.common.genome.refgenome;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class GeneNameMappingTest {

    private final GeneNameMapping victim = new GeneNameMapping();

    @Test
    public void testPRAMEF11() {
        assertTrue(victim.isValidV37Gene("PRAMEF11"));
        assertTrue(victim.isValidV37Gene("WI2-3308P17.2"));
        assertTrue(victim.isValidV38Gene("PRAMEF11"));

        assertFalse(victim.isValidV38Gene("WI2-3308P17.2"));
        assertEquals("WI2-3308P17.2", victim.v37Gene("PRAMEF11"));
    }

    @Test
    public void testPERM1() {
        assertTrue(victim.isValidV37Gene("C1orf170"));
        assertFalse(victim.isValidV37Gene("PERM1"));

        assertTrue(victim.isValidV38Gene("PERM1"));
        assertFalse(victim.isValidV38Gene("C1orf170"));

        assertEquals("C1orf170", victim.v37Gene("PERM1"));
    }

    @Test
    public void testPOTE() {
        assertTrue(victim.isValidV37Gene("POTEM"));
        assertTrue(victim.isValidV37Gene("POTEG"));
        assertTrue(victim.isValidV37Gene("POTEH"));

        assertTrue(victim.isValidV38Gene("POTEM"));
        assertTrue(victim.isValidV38Gene("POTEG"));
        assertTrue(victim.isValidV38Gene("POTEH"));

        assertEquals("POTEM", victim.v37Gene("POTEG"));
        assertEquals("POTEG", victim.v37Gene("POTEM"));
        assertEquals("POTEH", victim.v37Gene("POTEH"));
    }
}

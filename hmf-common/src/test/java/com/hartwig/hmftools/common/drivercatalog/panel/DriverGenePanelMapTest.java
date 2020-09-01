package com.hartwig.hmftools.common.drivercatalog.panel;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class DriverGenePanelMapTest {

    DriverGenePanelMap victim = new DriverGenePanelMap();

    @Test
    public void testPRAMEF11() {
        assertTrue(victim.isValidHg19Gene("PRAMEF11"));
        assertTrue(victim.isValidHg19Gene("WI2-3308P17.2"));
        assertTrue(victim.isValidHg38Gene("PRAMEF11"));

        assertFalse(victim.isValidHg38Gene("WI2-3308P17.2"));
        assertEquals("WI2-3308P17.2", victim.hg19Gene("PRAMEF11"));

    }

    @Test
    public void testPERM1() {
        assertTrue(victim.isValidHg19Gene("C1orf170"));
        assertFalse(victim.isValidHg19Gene("PERM1"));

        assertTrue(victim.isValidHg38Gene("PERM1"));
        assertFalse(victim.isValidHg38Gene("C1orf170"));

        assertEquals("C1orf170", victim.hg19Gene("PERM1"));
    }

    @Test
    public void testPOTE() {
        assertTrue(victim.isValidHg19Gene("POTEM"));
        assertTrue(victim.isValidHg19Gene("POTEG"));
        assertTrue(victim.isValidHg19Gene("POTEH"));

        assertTrue(victim.isValidHg38Gene("POTEM"));
        assertTrue(victim.isValidHg38Gene("POTEG"));
        assertTrue(victim.isValidHg38Gene("POTEH"));

        assertEquals("POTEM", victim.hg19Gene("POTEG"));
        assertEquals("POTEG", victim.hg19Gene("POTEM"));
        assertEquals("POTEH", victim.hg19Gene("POTEH"));
    }
}

package com.hartwig.hmftools.common.genome.refgenome;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.genome.genepanel.GeneNameMapping37to38;

import org.junit.Test;

public class GeneNameMapping37to38Test
{

    private final GeneNameMapping37to38 victim = GeneNameMapping37to38.loadFromEmbeddedResource();

    @Test
    public void mappingWorksForPRAMEF11() {
        assertTrue(victim.isValidV37Gene("PRAMEF11"));
        assertTrue(victim.isValidV38Gene("PRAMEF11"));

        assertTrue(victim.isValidV37Gene("WI2-3308P17.2"));
        assertFalse(victim.isValidV38Gene("WI2-3308P17.2"));
        assertEquals("WI2-3308P17.2", victim.v37Gene("PRAMEF11"));
    }

    @Test
    public void mappingWorksForPERM1() {
        assertTrue(victim.isValidV37Gene("C1orf170"));
        assertFalse(victim.isValidV37Gene("PERM1"));

        assertTrue(victim.isValidV38Gene("PERM1"));
        assertFalse(victim.isValidV38Gene("C1orf170"));

        assertEquals("C1orf170", victim.v37Gene("PERM1"));
    }

    @Test
    public void mappingWorksForPOTE() {
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

    @Test
    public void promiscuousIGGenesAreRetained() {
        for (String igGene : GeneNameMapping37to38.IG_GENES) {
            assertTrue(victim.isValidV37Gene(igGene));
            assertTrue(victim.isValidV38Gene(igGene));

            assertEquals(igGene, victim.v37Gene(igGene));
            assertEquals(igGene, victim.v38Gene(igGene));
        }
    }
}

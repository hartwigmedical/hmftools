package com.hartwig.hmftools.serve.sources.vicc.curation;

import static org.junit.Assert.*;

import org.junit.Test;

public class FusionCurationTest {

    @Test
    public void canCurateFusions() {
        assertEquals("CCDC6-RET", FusionCuration.curatedFusions("RET-CCDC6"));
        assertNotEquals("PDGFRB-CAPRIN1", FusionCuration.curatedFusions("GPIAP1-PDGFRB"));
    }

}
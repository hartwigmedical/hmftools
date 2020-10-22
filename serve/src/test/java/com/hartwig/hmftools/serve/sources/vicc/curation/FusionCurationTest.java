package com.hartwig.hmftools.serve.sources.vicc.curation;

import static org.junit.Assert.*;

import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;

import org.junit.Ignore;
import org.junit.Test;

public class FusionCurationTest {

    @Test
    @Ignore
    public void canCurateFusions() {
        Feature feature = ImmutableFeature.builder().geneSymbol("a").name("a").build();

        assertEquals("CCDC6-RET", FusionCuration.curatedFusions("RET-CCDC6", feature));
        assertNotEquals("PDGFRB-CAPRIN1", FusionCuration.curatedFusions("GPIAP1-PDGFRB", feature));
    }

}
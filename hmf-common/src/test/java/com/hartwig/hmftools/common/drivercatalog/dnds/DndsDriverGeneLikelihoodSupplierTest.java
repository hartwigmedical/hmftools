package com.hartwig.hmftools.common.drivercatalog.dnds;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import org.junit.Before;
import org.junit.Test;

public class DndsDriverGeneLikelihoodSupplierTest {

    private Map<String, DndsDriverGeneLikelihood> tsg;
    private Map<String, DndsDriverImpactLikelihood> onco;

    @Before
    public void setup() {
        onco = DndsDriverGeneLikelihoodSupplier.oncoLikelihood();
        tsg = DndsDriverGeneLikelihoodSupplier.tsgLikelihood();
    }

    @Test
    public void testReadOncoGenes() {
        final DndsDriverImpactLikelihood missense = onco.get("AKT1");

        assertEquals(0.412, missense.dndsLikelihood(), 0.001);
        assertEquals(0.002, missense.pDriver(), 0.001);
        assertEquals(1e-07, missense.pVariantNonDriverFactor(), 1e-7);
    }

    @Test
    public void testReadTSGGenes() {

        final DndsDriverGeneLikelihood gene = tsg.get("ACVR1B");
        final DndsDriverImpactLikelihood missense = gene.missense();

        assertEquals("ACVR1B", gene.gene());
        assertEquals(0.420, missense.dndsLikelihood(), 0.001);
        assertEquals(0.003, missense.pDriver(), 0.001);
        assertEquals(2e-07, missense.pVariantNonDriverFactor(), 1e-7);
    }

    @Test
    public void testBiallelic() {
        final DndsDriverGeneLikelihood apc = tsg.get("APC");
        assertTrue(apc.useBiallelic());
        assertNotEquals(apc.missense(), apc.missenseBiallelic());
        assertNotEquals(apc.missense(), apc.missenseNonBiallelic());
        assertNotEquals(apc.missenseBiallelic(), apc.missenseNonBiallelic());

        final DndsDriverGeneLikelihood sox9 = tsg.get("SOX9");
        assertFalse(sox9.useBiallelic());
        assertEquals(sox9.missense(), sox9.missenseBiallelic());
        assertEquals(sox9.missense(), sox9.missenseNonBiallelic());
    }
}

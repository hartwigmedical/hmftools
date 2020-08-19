package com.hartwig.hmftools.common.drivercatalog.dnds;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;
import java.util.stream.Collectors;

import org.junit.Before;
import org.junit.Test;

public class DndsDriverGeneLikelihoodSupplierTest {

    private Map<String, DndsDriverGeneLikelihood> tsg;
    private Map<String, DndsDriverImpactLikelihood> onco;

    @Before
    public void setup() {
        onco = DndsDriverGeneLikelihoodSupplier.oncoLikelihood()
                .values()
                .stream()
                .collect(Collectors.toMap(DndsDriverGeneLikelihood::gene, DndsDriverGeneLikelihood::missense));
        tsg = DndsDriverGeneLikelihoodSupplier.tsgLikelihood();
    }

    @Test
    public void testReadOncoGenes() {
        DndsDriverImpactLikelihood missense = onco.get("ABL1");

        assertEquals(0.0002922, missense.driversPerSample(), 1e-7);
        assertEquals(3.88883831596225e-07, missense.passengersPerMutation(), 1e-8);
    }

    @Test
    public void testReadTSGGenes() {
        DndsDriverGeneLikelihood gene = tsg.get("ACVR1B");
        DndsDriverImpactLikelihood missense = gene.missense();

        assertEquals("ACVR1B", gene.gene());
        assertEquals(0.003, missense.driversPerSample(), 0.001);
        assertEquals(2e-07, missense.passengersPerMutation(), 1e-7);
    }

    @Test
    public void testBiallelic() {
        DndsDriverGeneLikelihood apc = tsg.get("APC");
        assertTrue(apc.useBiallelic());
        assertNotEquals(apc.missense(), apc.missenseBiallelic());
        assertNotEquals(apc.missense(), apc.missenseNonBiallelic());
        assertNotEquals(apc.missenseBiallelic(), apc.missenseNonBiallelic());

        DndsDriverGeneLikelihood sox9 = tsg.get("SOX9");
        assertFalse(sox9.useBiallelic());
        assertEquals(sox9.missense(), sox9.missenseBiallelic());
        assertEquals(sox9.missense(), sox9.missenseNonBiallelic());
    }
}

package com.hartwig.hmftools.common.dnds;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import org.junit.Test;

public class DndsDriverGeneLikelihoodSupplierTest {

    @Test
    public void testReadOncoGenes() {
        final Map<String, DndsDriverGeneLikelihood> dndsLiklihoods = DndsDriverGeneLikelihoodSupplier.oncoLikelihood();

        // ABL1	0.159495158514378	0.00192322644362452	4.09846358113461e-07	0	0	0	0	0	0	0	0	0
        final DndsDriverGeneLikelihood gene = dndsLiklihoods.get("ABL1");
        final DndsDriverImpactLikelihood missense = gene.missense();

        assertEquals("ABL1", gene.gene());
        assertEquals(0.159, missense.dndsLikelihood(), 0.001);
        assertEquals(0.002, missense.pDriver(), 0.001);
        assertEquals(4e-07, missense.pVariantNonDriverFactor(), 1e-7);
    }

    @Test
    public void testReadTSGGenes() {
        final Map<String, DndsDriverGeneLikelihood> dndsLiklihoods = DndsDriverGeneLikelihoodSupplier.tsgLikelihood();

        //ACVR1B	0.169657430054376	0.00148141622916503	2.93196959039058e-07	0.813987487614234	0.00169228167903167	1.56384975335264e-08	0.867133919724453	0.00108166393312822	6.70221322865416e-09	0.0923460717793221	7.67950700867543e-05	1.89129931736336e-07
        final DndsDriverGeneLikelihood gene = dndsLiklihoods.get("ACVR1B");
        final DndsDriverImpactLikelihood missense = gene.missense();

        assertEquals("ACVR1B", gene.gene());
        assertEquals(0.170, missense.dndsLikelihood(), 0.001);
        assertEquals(0.001, missense.pDriver(), 0.001);
        assertEquals(2e-07, missense.pVariantNonDriverFactor(), 1e-7);

    }

}

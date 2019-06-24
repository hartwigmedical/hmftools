package com.hartwig.hmftools.common.dnds;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import org.junit.Test;

public class DndsDriverGeneLikelihoodSupplierTest {

    @Test
    public void testReadOncoGenes() {
        final Map<String, DndsDriverGeneLikelihood> dndsLiklihoods = DndsDriverGeneLikelihoodSupplier.oncoLikelihood();

        // AKT1	0.202114162676839	0.000802950703029438	1.24219636244432e-07	0	0	0	0	0	0	0	0	0
        final DndsDriverGeneLikelihood gene = dndsLiklihoods.get("AKT1");
        final DndsDriverImpactLikelihood missense = gene.missense();

        assertEquals("AKT1", gene.gene());
        assertEquals(0.202, missense.dndsLikelihood(), 0.001);
        assertEquals(8.0e-04, missense.pDriver(), 0.001);
        assertEquals(1e-07, missense.pVariantNonDriverFactor(), 1e-7);
    }

    @Test
    public void testReadTSGGenes() {
        final Map<String, DndsDriverGeneLikelihood> dndsLikelihoods = DndsDriverGeneLikelihoodSupplier.tsgLikelihood();

        //ACVR1B
        // 0.402200555288277	0.00353808661008416	2.0608133653813e-07
        // 0.832357869758091	0.0014171813900535	1.11855123684052e-08
        // 0.856306745506935	0.0009719713342871	6.39172135337441e-09	0.44198710815399	0.000501687977473315	1.67959675513019e-07
        final DndsDriverGeneLikelihood gene = dndsLikelihoods.get("ACVR1B");
        final DndsDriverImpactLikelihood missense = gene.missense();

        assertEquals("ACVR1B", gene.gene());
        assertEquals(0.402, missense.dndsLikelihood(), 0.001);
        assertEquals(0.003, missense.pDriver(), 0.001);
        assertEquals(2e-07, missense.pVariantNonDriverFactor(), 1e-7);
    }
}

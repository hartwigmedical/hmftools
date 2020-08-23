package com.hartwig.hmftools.common.drivercatalog.dnds;

import static org.junit.Assert.assertEquals;

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

}

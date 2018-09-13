package com.hartwig.hmftools.common.dnds;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Map;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.junit.Test;

public class DndsDriverLikelihoodSupplierTest {

    private static final double PRECISION = 0.01;

    @Test
    public void testReadFromResource() throws IOException {
        final Map<String, DndsDriverLikelihood> dndsLiklihoods = DndsDriverLikelihoodSupplier.oncoLikelihood();
//        assertEquals(1, dndsLiklihoods.get("TP53").indel(), PRECISION);
//        assertEquals(0.87, dndsLiklihoods.get("TP53").missense(), PRECISION);
//        assertEquals(0.98, dndsLiklihoods.get("TP53").nonsense(), PRECISION);
//        assertEquals(1, dndsLiklihoods.get("TP53").splice(), PRECISION);
    }

    @Test
    public void testStuff() {
        PoissonDistribution myPoisson =  new PoissonDistribution(0.3);
        System.out.println(myPoisson.cumulativeProbability(0));

    }

}

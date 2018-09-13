package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Map;

import com.hartwig.hmftools.common.dnds.DndsDriverLikelihood;
import com.hartwig.hmftools.common.dnds.DndsDriverLikelihoodSupplier;

import org.junit.Test;

public class OncoDriversTest {

    private static final double EPSILON = 0.0001;


    @Test
    public void testHIST2H3D() throws IOException {
        final Map<String, DndsDriverLikelihood> dnds = DndsDriverLikelihoodSupplier.oncoLikelihood();
        double value = OncoDrivers.missenseProbabilityDriverVariant(27742, dnds.get("HIST2H3D"));
        assertEquals(0.7042, value, EPSILON);
    }

    @Test
    public void testABL1() throws IOException {
        final Map<String, DndsDriverLikelihood> dnds = DndsDriverLikelihoodSupplier.oncoLikelihood();
        double value = OncoDrivers.missenseProbabilityDriverVariant(996698, dnds.get("ABL1"));
        assertEquals(0.0057, value, EPSILON);
    }
}

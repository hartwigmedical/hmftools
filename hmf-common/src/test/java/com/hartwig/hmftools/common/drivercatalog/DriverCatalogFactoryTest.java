package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;

import org.junit.Test;

public class DriverCatalogFactoryTest {

    private static final double EPSILON = 0.0001;

    @Test
    public void testHIST2H3DMissense() {
        final Map<String, DndsDriverGeneLikelihood> dnds = DndsDriverGeneLikelihoodSupplier.oncoLikelihood();
        double value = DriverCatalogFactory.probabilityDriverVariant(27742, dnds.get("HIST2H3D").missense());
        assertEquals(0.7042, value, EPSILON);
    }

    @Test
    public void testABL1Missense()  {
        final Map<String, DndsDriverGeneLikelihood> dnds = DndsDriverGeneLikelihoodSupplier.oncoLikelihood();
        double value = DriverCatalogFactory.probabilityDriverVariant(996698, dnds.get("ABL1").missense());
        assertEquals(0.0057, value, EPSILON);
    }

    @Test
    public void testGATA3Indel()  {
        final Map<String, DndsDriverGeneLikelihood> dnds = DndsDriverGeneLikelihoodSupplier.tsgLikelihood();
        double value = DriverCatalogFactory.probabilityDriverVariant(587, dnds.get("GATA3").indel());
        assertEquals(0.9952, value, EPSILON);
    }
}

package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.dnds.ImmutableDndsDriverImpactLikelihood;

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

    @Test
    public void testMultipleZeroNonDriver()  {
        DndsDriverImpactLikelihood indelLikelihood = ImmutableDndsDriverImpactLikelihood.builder().dndsLikelihood(1).pDriver(0.01).pVariantNonDriverFactor(0).build();
        double value = DriverCatalogFactory.probabilityDriverVariant(1000, 1000, indelLikelihood, indelLikelihood);
        assertEquals(1, value, EPSILON);
    }

    @Test
    public void testZeroNonDriverWithStandard()  {
        DndsDriverImpactLikelihood indelLikelihood = ImmutableDndsDriverImpactLikelihood.builder().dndsLikelihood(1).pDriver(0.01).pVariantNonDriverFactor(0).build();
        DndsDriverImpactLikelihood missenseLikelihood = ImmutableDndsDriverImpactLikelihood.builder().dndsLikelihood(0.87).pDriver(0.01).pVariantNonDriverFactor(10e-8).build();

        double value = DriverCatalogFactory.probabilityDriverVariant(10000, 1000, missenseLikelihood, indelLikelihood);
        assertEquals(1, value, EPSILON);
    }

}

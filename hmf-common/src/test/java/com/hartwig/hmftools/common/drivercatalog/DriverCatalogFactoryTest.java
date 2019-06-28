package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.dnds.ImmutableDndsDriverImpactLikelihood;

import org.junit.Before;
import org.junit.Test;

public class DriverCatalogFactoryTest {

    private static final double EPSILON = 0.0001;

    private Map<String, DndsDriverGeneLikelihood> tsg;
    private Map<String, DndsDriverImpactLikelihood> onco;

    @Before
    public void setup() {
        onco = DndsDriverGeneLikelihoodSupplier.oncoLikelihood();
        tsg = DndsDriverGeneLikelihoodSupplier.tsgLikelihood();
    }

    @Test
    public void testHIST2H3DMissense() {
        double value = DriverCatalogFactory.probabilityDriverVariant(27742, onco.get("HIST2H3D"));
        assertEquals(0.6065, value, EPSILON);
    }

    @Test
    public void testABL1Missense()  {
        double value = DriverCatalogFactory.probabilityDriverVariant(996698, onco.get("ABL1"));
        assertEquals(0.0009, value, EPSILON);
    }

    @Test
    public void testGATA3Indel()  {
        double value = DriverCatalogFactory.probabilityDriverVariant(587, tsg.get("GATA3").indel());
        assertEquals(0.9945, value, EPSILON);
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

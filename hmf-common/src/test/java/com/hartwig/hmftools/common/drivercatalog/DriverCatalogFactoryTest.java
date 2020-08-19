package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.ImmutableDndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelFactory;

import org.junit.Before;
import org.junit.Test;

public class DriverCatalogFactoryTest {

    private static final double EPSILON = 0.0001;

    private Map<String, DndsDriverGeneLikelihood> tsg;
    private Map<String, DndsDriverImpactLikelihood> onco;

    @Before
    public void setup() {
        DriverGenePanel genePanel = new DriverGenePanelFactory().create();
        onco = genePanel.oncoLikelihood()
                .values()
                .stream()
                .collect(Collectors.toMap(DndsDriverGeneLikelihood::gene, DndsDriverGeneLikelihood::missense));
        tsg = genePanel.tsgLikelihood();
    }

    @Test
    public void testHIST2H3DMissense() {
        double value = DriverCatalogFactory.probabilityDriverVariant(27742, onco.get("HIST2H3D"));
        assertEquals(0.6065, value, EPSILON);
    }

    @Test
    public void testABL1Missense() {
        double value = DriverCatalogFactory.probabilityDriverVariant(996698, onco.get("ABL1"));
        assertEquals(0.0009, value, EPSILON);
    }

    @Test
    public void testGATA3Indel() {
        double value = DriverCatalogFactory.probabilityDriverVariant(587, tsg.get("GATA3").indel());
        assertEquals(0.9945, value, EPSILON);
    }

    @Test
    public void testMultipleZeroNonDriver() {
        DndsDriverImpactLikelihood indelLikelihood =
                ImmutableDndsDriverImpactLikelihood.builder().dndsLikelihood(1).pDriver(0.01).pVariantNonDriverFactor(0).build();
        double value = DriverCatalogFactory.probabilityDriverVariant(1000, 1000, indelLikelihood, indelLikelihood);
        assertEquals(0, value, EPSILON);
    }

    @Test
    public void testFallBackOnSingleProbabilityIfMultiFailsDueToZeroValues() {
        DndsDriverImpactLikelihood nonsense =
                ImmutableDndsDriverImpactLikelihood.builder().dndsLikelihood(0.7).pDriver(2E-4).pVariantNonDriverFactor(4E-9).build();
        double singleNonsenseLikelihood = DriverCatalogFactory.probabilityDriverVariant(1, nonsense);
        assertTrue(singleNonsenseLikelihood - 0.1 > nonsense.dndsLikelihood());

        DndsDriverImpactLikelihood splice =
                ImmutableDndsDriverImpactLikelihood.builder().dndsLikelihood(0).pDriver(0).pVariantNonDriverFactor(0).build();
        assertEquals(0, splice.dndsLikelihood(), EPSILON);
        assertEquals(0, DriverCatalogFactory.probabilityDriverVariant(1, splice), EPSILON);

        double victim = DriverCatalogFactory.probabilityDriverVariant(1, 1, nonsense, splice);
        assertEquals(singleNonsenseLikelihood, victim, EPSILON);

        victim = DriverCatalogFactory.probabilityDriverVariant(1, 1, splice, nonsense);
        assertEquals(singleNonsenseLikelihood, victim, EPSILON);
    }

    @Test
    public void testZeroNonDriverWithStandard() {
        DndsDriverImpactLikelihood indelLikelihood =
                ImmutableDndsDriverImpactLikelihood.builder().dndsLikelihood(1).pDriver(0.01).pVariantNonDriverFactor(0).build();
        DndsDriverImpactLikelihood missenseLikelihood =
                ImmutableDndsDriverImpactLikelihood.builder().dndsLikelihood(0.87).pDriver(0.01).pVariantNonDriverFactor(10e-8).build();

        double expectedMissense = DriverCatalogFactory.probabilityDriverVariant(10000, missenseLikelihood);
        double value = DriverCatalogFactory.probabilityDriverVariant(10000, 1000, missenseLikelihood, indelLikelihood);
        assertEquals(expectedMissense, value, EPSILON);
    }
}

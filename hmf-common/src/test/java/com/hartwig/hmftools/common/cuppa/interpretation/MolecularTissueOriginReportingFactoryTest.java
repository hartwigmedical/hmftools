package com.hartwig.hmftools.common.cuppa.interpretation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import org.junit.Test;

public class MolecularTissueOriginReportingFactoryTest  {

    private static final double EPSILON = 1.0E-10;

    @Test
    public void canCreateCuppaReportingKnown() {
        CuppaPrediction prediction = ImmutableCuppaPrediction.builder().cancerType("Melanoma").likelihood(0.90).build();
        MolecularTissueOriginReporting reporting = MolecularTissueOriginReportingFactory.createMolecularTissueOriginReportingData(prediction);
        assertEquals("Melanoma", reporting.bestCancerType());
        assertEquals(0.90, reporting.bestLikelihood(), EPSILON);
        assertEquals("Melanoma", reporting.interpretCancerType());
        assertEquals(0.90, reporting.interpretLikelihood(), EPSILON);
    }

    @Test
    public void canCreateCuppaReportingKnownUterus() {
        CuppaPrediction prediction = ImmutableCuppaPrediction.builder().cancerType("Uterus: Endometrium").likelihood(0.90).build();
        MolecularTissueOriginReporting reporting = MolecularTissueOriginReportingFactory.createMolecularTissueOriginReportingData(prediction);
        assertEquals("Endometrium", reporting.bestCancerType());
        assertEquals(0.90, reporting.bestLikelihood(), EPSILON);
        assertEquals("Endometrium", reporting.interpretCancerType());
        assertEquals(0.90, reporting.interpretLikelihood(), EPSILON);
    }

    @Test
    public void canCreateCuppaReportingKnownColon() {
        CuppaPrediction prediction = ImmutableCuppaPrediction.builder().cancerType("Colorectum/Appendix/SmallIntestine").likelihood(0.90).build();
        MolecularTissueOriginReporting reporting = MolecularTissueOriginReportingFactory.createMolecularTissueOriginReportingData(prediction);
        assertEquals("Lower GI tract", reporting.bestCancerType());
        assertEquals(0.90, reporting.bestLikelihood(), EPSILON);
        assertEquals("Lower GI tract", reporting.interpretCancerType());
        assertEquals(0.90, reporting.interpretLikelihood(), EPSILON);
    }

    @Test
    public void canCreateCuppaReportingInconclusiveWithLikelihood() {
        CuppaPrediction prediction = ImmutableCuppaPrediction.builder().cancerType("Melanoma").likelihood(0.60).build();
        MolecularTissueOriginReporting reporting = MolecularTissueOriginReportingFactory.createMolecularTissueOriginReportingData(prediction);
        assertEquals("Melanoma", reporting.bestCancerType());
        assertEquals(0.60, reporting.bestLikelihood(), EPSILON);
        assertEquals("results inconclusive", reporting.interpretCancerType());
        assertNull(reporting.interpretLikelihood());
    }

    @Test
    public void canCreateCuppaReportingInconclusiveWithoutLikelihood() {
        CuppaPrediction prediction = ImmutableCuppaPrediction.builder().cancerType("Melanoma").likelihood(0.40).build();
        MolecularTissueOriginReporting reporting = MolecularTissueOriginReportingFactory.createMolecularTissueOriginReportingData(prediction);
        assertEquals("Melanoma", reporting.bestCancerType());
        assertEquals(0.40, reporting.bestLikelihood(), EPSILON);
        assertEquals("results inconclusive", reporting.interpretCancerType());
        assertNull(reporting.interpretLikelihood());
    }
}
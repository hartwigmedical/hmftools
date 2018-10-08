package com.hartwig.hmftools.patientreporter.algo;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.hartwig.hmftools.patientreporter.PatientReporterTestUtil;

import org.junit.Test;

public class DrupActionabilityModelFactoryTest {

    @Test
    public void canReadTestDrupCsv() throws IOException {
        DrupActionabilityModel drupActionabilityModel = PatientReporterTestUtil.testDrupActionabilityModel();

        assertEquals(2, drupActionabilityModel.actionableGenes().size());
        assertEquals(1, drupActionabilityModel.geneDriverCategoryMap().size());
    }
}
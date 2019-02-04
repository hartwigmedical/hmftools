package com.hartwig.hmftools.patientreporter.structural;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class FusionFactoryTest {

    private static final String FUSION_FILE = Resources.getResource("test_run/svAnalysis/CPCT11111111T_fusions.csv").getPath();

    @Test
    public void canReadFromFile() throws IOException {
       FusionAnalyzer analyzer = FusionFactory.fromFusionFile(FUSION_FILE);
       assertEquals(1, analyzer.fusions().size());
    }
}

package com.hartwig.hmftools.patientreporter.pipeline;

import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReporter;

import org.junit.Test;

public class PipelineVersionTest  {

    @Test
    public void testPipelineVersion() {
        PipelineVersion.checkPipelineVersion("5.22", "5.22", false);
        PipelineVersion.checkPipelineVersion("5.22", "5.22", true);
        PipelineVersion.checkPipelineVersion("5.22", "5.21", true);
        PipelineVersion.checkPipelineVersion(null, "5.21", true);
    }

    @Test(expected = IllegalArgumentException.class)
    public void crashTestPipelineVersionOnNoActual() {
        PipelineVersion.checkPipelineVersion(null, "5.21", false);
    }

    @Test(expected = IllegalArgumentException.class)
    public void crashTestPipelineVersionOnVersionDifference() {
        PipelineVersion.checkPipelineVersion("5.22", "5.21", false);
    }

}
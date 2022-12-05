package com.hartwig.hmftools.patientdb.clinical.lims;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class LimsAnalysisTypeTest {

    @Test
    public void canExtractLimsAnalysisType() {
        assertEquals(LimsAnalysisType.SOMATIC_R, LimsAnalysisType.extractAnalysisType("Somatic_R"));
        assertEquals(LimsAnalysisType.SOMATIC_T, LimsAnalysisType.extractAnalysisType("Somatic_T"));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownAnalysisType() {
        //noinspection ResultOfMethodCallIgnored
        LimsAnalysisType.extractAnalysisType("somatic_T");
    }
}
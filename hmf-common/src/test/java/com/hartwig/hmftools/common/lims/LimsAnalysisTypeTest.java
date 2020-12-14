package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.*;

import org.junit.Test;

public class LimsAnalysisTypeTest {
    @Test
    public void canExtractLimsAnalysisType() {
        assertEquals(LimsAnalysisType.SOMATIC_R, LimsAnalysisType.extractAnalysisType("Somatic_R"));
        assertEquals(LimsAnalysisType.SOMATIC_T, LimsAnalysisType.extractAnalysisType("Somatic_T"));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownCohortType() {
        LimsAnalysisType.extractAnalysisType("somatic_T");
    }

}
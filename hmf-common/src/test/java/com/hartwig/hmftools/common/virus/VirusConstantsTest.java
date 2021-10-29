package com.hartwig.hmftools.common.virus;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.apache.commons.compress.utils.Lists;
import org.junit.Test;

public class VirusConstantsTest {

    @Test
    public void canExtractVirusConstants() {
        VirusConstants VirusConstantsHPV = VirusConstants.fromVirusName("HPV");
        assertEquals(VirusConstants.HPV, VirusConstantsHPV);
        assertTrue(VirusConstantsHPV.reportVirusOnSummary());

        VirusConstants VirusConstantsMCV = VirusConstants.fromVirusName("MCV");
        assertEquals(VirusConstants.MCV, VirusConstantsMCV);
        assertTrue(VirusConstantsMCV.reportVirusOnSummary());

        VirusConstants VirusConstantsEBV = VirusConstants.fromVirusName("EBV");
        assertEquals(VirusConstants.EBV, VirusConstantsEBV);
        assertTrue(VirusConstantsEBV.reportVirusOnSummary());

        VirusConstants VirusConstantsHBV = VirusConstants.fromVirusName("HBV");
        assertEquals(VirusConstants.HBV, VirusConstantsHBV);
        assertTrue(VirusConstantsHBV.reportVirusOnSummary());

        VirusConstants VirusConstantsHHV8 = VirusConstants.fromVirusName("HHV-8");
        assertEquals(VirusConstants.HHV8, VirusConstantsHHV8);
        assertTrue(VirusConstantsHHV8.reportVirusOnSummary());
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownVirusConstants() {
        //noinspection ResultOfMethodCallIgnored
        VirusConstants.fromVirusName("ABC");
    }

    @Test
    public void canExtractAllSummaryReportedViruses() {
        List<String> reportableSummaryViruses = VirusConstants.allViruses();
        List<String> expectedViruses = Lists.newArrayList();
        expectedViruses.add("MCV");
        expectedViruses.add("EBV");
        expectedViruses.add("HPV");
        expectedViruses.add("HBV");
        expectedViruses.add("HHV-8");
        assertEquals(reportableSummaryViruses, expectedViruses);
    }
}
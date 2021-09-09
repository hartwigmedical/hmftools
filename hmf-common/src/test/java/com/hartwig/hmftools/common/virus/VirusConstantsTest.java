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
        VirusConstants VirusConstantsHPV = VirusConstants.virusName("HPV");
        assertEquals(VirusConstants.HPV, VirusConstantsHPV);
        assertTrue(VirusConstants.virusName("HPV").reportedVirusOnSummary());

        VirusConstants VirusConstantsMCV = VirusConstants.virusName("MCV");
        assertEquals(VirusConstants.MCV, VirusConstantsMCV);
        assertTrue(VirusConstants.virusName("MCV").reportedVirusOnSummary());

        VirusConstants VirusConstantsEBV = VirusConstants.virusName("EBV");
        assertEquals(VirusConstants.EBV, VirusConstantsEBV);
        assertTrue(VirusConstants.virusName("EBV").reportedVirusOnSummary());

        VirusConstants VirusConstantsHBV = VirusConstants.virusName("HBV");
        assertEquals(VirusConstants.HBV, VirusConstantsHBV);
        assertFalse(VirusConstants.virusName("HBV").reportedVirusOnSummary());

        VirusConstants VirusConstantsHHV8 = VirusConstants.virusName("HHV-8");
        assertEquals(VirusConstants.HHV8, VirusConstantsHHV8);
        assertFalse(VirusConstants.virusName("HHV-8").reportedVirusOnSummary());
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownVirusConstants() {
        VirusConstants.virusName("ABC");
    }

    @Test
    public void canExtractAllSummaryReportedVirussen() {
        List<String> reportableSummaryVirussen = VirusConstants.allVirussen();
        List<String> expectedVirussen = Lists.newArrayList();
        expectedVirussen.add("MCV");
        expectedVirussen.add("EBV");
        expectedVirussen.add("HPV");
        assertEquals(reportableSummaryVirussen, expectedVirussen);
    }
}
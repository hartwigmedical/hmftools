package com.hartwig.hmftools.patientreporter.report.util;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class PatientReportFormatTest {

    @Test
    public void removesCopiesWhenNoTumor() {
        assertTrue(PatientReportFormat.correctValueForFitReliability("1", true).equals("1"));
        assertTrue(PatientReportFormat.correctValueForFitReliability("1", false).equals("N/A"));
    }

    @Test
    public void negativeMinorAllelePloidyWorks() {
        assertEquals("AAAAAAA", PatientReportFormat.descriptiveBAF(7, -1));
    }

    @Test
    public void descriptiveBAFWorksAroundBoundary() {
        assertEquals("AA", PatientReportFormat.descriptiveBAF(2, 2));
        assertEquals("", PatientReportFormat.descriptiveBAF(-12, 0));
    }

    @Test
    public void descriptiveBAFWorksForTypicalCases() {
        assertEquals("AA", PatientReportFormat.descriptiveBAF(2, 0));
        assertEquals("AB", PatientReportFormat.descriptiveBAF(2, 1));
        assertEquals("AB", PatientReportFormat.descriptiveBAF(2, 1.4));
        assertEquals("AA", PatientReportFormat.descriptiveBAF(2, 1.6));
        assertEquals("AA", PatientReportFormat.descriptiveBAF(2, 2));

        assertEquals("AAA", PatientReportFormat.descriptiveBAF(3, 0));
        assertEquals("AAB", PatientReportFormat.descriptiveBAF(3, 1.5));
        assertEquals("AAB", PatientReportFormat.descriptiveBAF(3, 2.49));
        assertEquals("AAA", PatientReportFormat.descriptiveBAF(3, 2.51));
        assertEquals("AAA", PatientReportFormat.descriptiveBAF(3, 3));

        assertEquals("AAAA", PatientReportFormat.descriptiveBAF(4, 0));
        assertEquals("AAAA", PatientReportFormat.descriptiveBAF(4, 0.1));
        assertEquals("AABB", PatientReportFormat.descriptiveBAF(4, 2));
        assertEquals("AABB", PatientReportFormat.descriptiveBAF(4, 2.49));
        assertEquals("AAAB", PatientReportFormat.descriptiveBAF(4, 2.51));
        assertEquals("AAAB", PatientReportFormat.descriptiveBAF(4, 3.49));
        assertEquals("AAAA", PatientReportFormat.descriptiveBAF(4, 3.51));
        assertEquals("AAAA", PatientReportFormat.descriptiveBAF(4, 3.8));
        assertEquals("AAAA", PatientReportFormat.descriptiveBAF(4, 4));

        assertEquals("A[12x]", PatientReportFormat.descriptiveBAF(12, 0));
    }
}
package com.hartwig.hmftools.viridian.reporting;

import static com.hartwig.hmftools.viridian.reporting.VirusReportStatus.CLONALITY_UNEVALUABLE;
import static com.hartwig.hmftools.viridian.reporting.VirusReportStatus.NOT_CLONAL;
import static com.hartwig.hmftools.viridian.reporting.VirusReportStatus.NOT_REPORTABLE;
import static com.hartwig.hmftools.viridian.reporting.VirusReportStatus.REPORTED_ON_CLONALITY;
import static com.hartwig.hmftools.viridian.reporting.VirusReportStatus.REPORTED_ON_INTEGRATION;
import static com.hartwig.hmftools.viridian.reporting.VirusReportStatus.REPORTED_ON_INTEGRATION_AND_CLONALITY;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.viridian.detection.DetectedVirus;
import com.hartwig.hmftools.viridian.detection.common.GenomeStats;
import com.hartwig.hmftools.viridian.detection.common.SummaryStats;
import com.hartwig.hmftools.viridian.detection.select.RepresentativeStatus;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.OncologyGroupInfo;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.junit.Test;

public class VirusReporterTest
{
    private static final OncologyGroup GROUP = new OncologyGroup("HPV16");
    private static final VirusGenome GENOME = new VirusGenome("c1", 1000, "HPV type 16", GROUP);
    private static final OncologyGroupInfo REPORTABLE_INFO = new OncologyGroupInfo(GROUP, VirusType.HPV, VirusLikelihoodType.HIGH);
    private static final OncologyGroupInfo NOT_REPORTABLE_INFO = new OncologyGroupInfo(GROUP, null, null);
    private static final double DEPTH = 10.0;
    private static final Double PER_COPY = 20.0;

    private static final boolean REPORTABLE = true;

    // The whole decision table. Thresholds are I=1 integration, C=0.5 copies per tumor cell. A null
    // copies-per-cell means clonality could not be evaluated.

    @Test
    public void testNotReportableAlwaysLosesRegardlessOfEvidence()
    {
        assertEquals(NOT_REPORTABLE, VirusReporter.decideStatus(!REPORTABLE, 5, 10.0));
        assertEquals(NOT_REPORTABLE, VirusReporter.decideStatus(!REPORTABLE, 0, null));
    }

    @Test
    public void testIntegrationReports()
    {
        assertEquals(REPORTED_ON_INTEGRATION, VirusReporter.decideStatus(REPORTABLE, 1, null));
    }

    @Test
    public void testIntegrationWithoutClonalLoad()
    {
        assertEquals(REPORTED_ON_INTEGRATION, VirusReporter.decideStatus(REPORTABLE, 1, 0.0));
    }

    @Test
    public void testIntegrationAndClonalityReportsAsBoth()
    {
        assertEquals(REPORTED_ON_INTEGRATION_AND_CLONALITY, VirusReporter.decideStatus(REPORTABLE, 1, 0.5));
    }

    @Test
    public void testClonalityReportsAtAndAboveFloor()
    {
        assertEquals(REPORTED_ON_CLONALITY, VirusReporter.decideStatus(REPORTABLE, 0, 0.5));
        assertEquals(REPORTED_ON_CLONALITY, VirusReporter.decideStatus(REPORTABLE, 0, 2.0));
    }

    @Test
    public void testClonalityBelowFloor()
    {
        assertEquals(NOT_CLONAL, VirusReporter.decideStatus(REPORTABLE, 0, 0.49));
    }

    @Test
    public void testClonalityUnevaluable()
    {
        assertEquals(CLONALITY_UNEVALUABLE, VirusReporter.decideStatus(REPORTABLE, 0, null));
    }

    @Test
    public void testReportedFlags()
    {
        assertTrue(REPORTED_ON_INTEGRATION.isReported());
        assertTrue(REPORTED_ON_CLONALITY.isReported());
        assertTrue(REPORTED_ON_INTEGRATION_AND_CLONALITY.isReported());
        assertFalse(NOT_REPORTABLE.isReported());
        assertFalse(NOT_CLONAL.isReported());
        assertFalse(CLONALITY_UNEVALUABLE.isReported());
    }

    // Assembling one group's report from detection, integration count, policy and the per-copy depth.

    // Representative depth 10 against per-copy depth 20 is 0.5 copies per tumor cell, which clears C.
    @Test
    public void testDetectedReportsOnClonality()
    {
        VirusReport report = VirusReporter.reportGroup(GROUP, detected(DEPTH), 0, REPORTABLE_INFO, PER_COPY);
        assertTrue(report.isDetected());
        assertEquals(0.5, report.copiesPerTumorCell(), 1e-9);
        assertEquals(REPORTED_ON_CLONALITY, report.status());
        assertTrue(report.isReported());
    }

    // Clonality needs the per-copy depth, which an unusable Purple fit leaves absent.
    @Test
    public void testDetectedWithoutPurpleIsUnevaluable()
    {
        VirusReport report = VirusReporter.reportGroup(GROUP, detected(DEPTH), 0, REPORTABLE_INFO, null);
        assertNull(report.copiesPerTumorCell());
        assertEquals(CLONALITY_UNEVALUABLE, report.status());
    }

    // An integration reports the virus even when it was not detected.
    @Test
    public void testIntegratedButNotDetected()
    {
        VirusReport report = VirusReporter.reportGroup(GROUP, null, 1, REPORTABLE_INFO, PER_COPY);
        assertFalse(report.isDetected());
        assertNull(report.copiesPerTumorCell());
        assertEquals(REPORTED_ON_INTEGRATION, report.status());
    }

    @Test
    public void testNotReportableIsNotReported()
    {
        VirusReport report = VirusReporter.reportGroup(GROUP, detected(DEPTH), 5, NOT_REPORTABLE_INFO, PER_COPY);
        assertEquals(NOT_REPORTABLE, report.status());
        assertFalse(report.isReported());
    }

    private static DetectedVirus detected(double meanDepth)
    {
        SummaryStats depth = SummaryStats.from(new int[] { (int) meanDepth });
        SummaryStats alignerScore = SummaryStats.from(new int[] { 30 });
        GenomeStats stats = new GenomeStats(GENOME, 1, 0, GENOME.length(), depth, alignerScore);
        return new DetectedVirus(stats, RepresentativeStatus.PASS);
    }
}

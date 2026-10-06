package com.hartwig.hmftools.viridian.reporting;

import static com.hartwig.hmftools.viridian.detection.select.OncologyGroupOutcome.MUTUAL;
import static com.hartwig.hmftools.viridian.detection.select.OncologyGroupOutcome.RESOLVED_CANDIDATES;
import static com.hartwig.hmftools.viridian.detection.select.OncologyGroupResolution.RESOLVED;
import static com.hartwig.hmftools.viridian.detection.select.OncologyGroupResolution.UNRESOLVED;
import static com.hartwig.hmftools.viridian.reporting.VirusReportStatus.COPY_NUMBER_TOO_LOW;
import static com.hartwig.hmftools.viridian.reporting.VirusReportStatus.COPY_NUMBER_UNEVALUABLE;
import static com.hartwig.hmftools.viridian.reporting.VirusReportStatus.NOT_REPORTABLE;
import static com.hartwig.hmftools.viridian.reporting.VirusReportStatus.REPORTED_ON_COPY_NUMBER;
import static com.hartwig.hmftools.viridian.reporting.VirusReportStatus.REPORTED_ON_INTEGRATION;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.viridian.detection.DetectedVirus;
import com.hartwig.hmftools.viridian.detection.common.ContigStats;
import com.hartwig.hmftools.viridian.detection.common.SummaryStats;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.OncologyGroupInfo;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import org.junit.Test;

public class VirusReporterTest
{
    private static final OncologyGroup GROUP = new OncologyGroup("HPV16");
    private static final ViralContig CONTIG = new ViralContig("c1", 1000, "HPV type 16", GROUP);
    private static final OncologyGroupInfo REPORTABLE_INFO = new OncologyGroupInfo(GROUP, VirusType.HPV, VirusLikelihoodType.HIGH);
    private static final OncologyGroupInfo NOT_REPORTABLE_INFO = new OncologyGroupInfo(GROUP, null, null);
    private static final double DEPTH = 10.0;
    private static final Double PER_COPY = 20.0;

    private static final boolean REPORTABLE = true;

    // The whole decision table. Thresholds are I=1 integration, C=0.5 copies per tumor cell. A null
    // copies-per-cell means copy number could not be evaluated.

    @Test
    public void testNotReportableAlwaysLosesRegardlessOfEvidence()
    {
        assertEquals(NOT_REPORTABLE, VirusReporter.decide(!REPORTABLE, 5, 10.0));
        assertEquals(NOT_REPORTABLE, VirusReporter.decide(!REPORTABLE, 0, null));
    }

    @Test
    public void testIntegrationReports()
    {
        assertEquals(REPORTED_ON_INTEGRATION, VirusReporter.decide(REPORTABLE, 1, null));
    }

    // Integration path is checked before copy number, so it wins even when copy number is too low.
    @Test
    public void testIntegrationTakesPrecedenceOverCopyNumber()
    {
        assertEquals(REPORTED_ON_INTEGRATION, VirusReporter.decide(REPORTABLE, 1, 0.0));
    }

    @Test
    public void testCopyNumberReportsAtAndAboveFloor()
    {
        assertEquals(REPORTED_ON_COPY_NUMBER, VirusReporter.decide(REPORTABLE, 0, 0.5));
        assertEquals(REPORTED_ON_COPY_NUMBER, VirusReporter.decide(REPORTABLE, 0, 2.0));
    }

    @Test
    public void testCopyNumberBelowFloor()
    {
        assertEquals(COPY_NUMBER_TOO_LOW, VirusReporter.decide(REPORTABLE, 0, 0.49));
    }

    @Test
    public void testCopyNumberUnevaluable()
    {
        assertEquals(COPY_NUMBER_UNEVALUABLE, VirusReporter.decide(REPORTABLE, 0, null));
    }

    @Test
    public void testReportedFlags()
    {
        assertTrue(REPORTED_ON_INTEGRATION.isReported());
        assertTrue(REPORTED_ON_COPY_NUMBER.isReported());
        assertFalse(NOT_REPORTABLE.isReported());
        assertFalse(COPY_NUMBER_TOO_LOW.isReported());
        assertFalse(COPY_NUMBER_UNEVALUABLE.isReported());
    }

    // Assembling one group's report from detection, integration count, policy and the per-copy depth.

    // Representative depth 10 against per-copy depth 20 is 0.5 copies per tumor cell, which clears C.
    @Test
    public void testPresentResolvedReportsOnCopyNumber()
    {
        VirusReport report = VirusReporter.reportGroup(GROUP, resolved(DEPTH), 0, REPORTABLE_INFO, PER_COPY);
        assertTrue(report.isPresent());
        assertEquals(0.5, report.copiesPerTumorCell(), 1e-9);
        assertEquals(REPORTED_ON_COPY_NUMBER, report.reason());
        assertTrue(report.isReported());
    }

    // Copy number needs the per-copy depth, which an unusable Purple fit leaves absent.
    @Test
    public void testPresentResolvedWithoutPurpleIsUnevaluable()
    {
        VirusReport report = VirusReporter.reportGroup(GROUP, resolved(DEPTH), 0, REPORTABLE_INFO, null);
        assertNull(report.copiesPerTumorCell());
        assertEquals(COPY_NUMBER_UNEVALUABLE, report.reason());
    }

    // Present but unresolved: no representative to measure, so copy number is unevaluable.
    @Test
    public void testPresentUnresolvedIsUnevaluable()
    {
        DetectedVirus unresolved = new DetectedVirus(GROUP, UNRESOLVED, MUTUAL, 100, 2, 2, 2, null);
        VirusReport report = VirusReporter.reportGroup(GROUP, unresolved, 0, REPORTABLE_INFO, PER_COPY);
        assertTrue(report.isPresent());
        assertNull(report.copiesPerTumorCell());
        assertEquals(COPY_NUMBER_UNEVALUABLE, report.reason());
    }

    // An integration reports the virus even when it was not called present.
    @Test
    public void testIntegratedButNotPresent()
    {
        VirusReport report = VirusReporter.reportGroup(GROUP, null, 1, REPORTABLE_INFO, PER_COPY);
        assertFalse(report.isPresent());
        assertNull(report.copiesPerTumorCell());
        assertEquals(REPORTED_ON_INTEGRATION, report.reason());
    }

    @Test
    public void testNotReportableIsNotReported()
    {
        VirusReport report = VirusReporter.reportGroup(GROUP, resolved(DEPTH), 5, NOT_REPORTABLE_INFO, PER_COPY);
        assertEquals(NOT_REPORTABLE, report.reason());
        assertFalse(report.isReported());
    }

    private static DetectedVirus resolved(double meanDepth)
    {
        SummaryStats depth = SummaryStats.from(new int[] { (int) meanDepth });
        SummaryStats alignerScore = SummaryStats.from(new int[] { 30 });
        ContigStats stats = new ContigStats(CONTIG, 1, 0, CONTIG.length(), depth, alignerScore);
        return new DetectedVirus(GROUP, RESOLVED, RESOLVED_CANDIDATES, 100, 1, 1, 0, stats);
    }
}

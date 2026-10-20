package com.hartwig.hmftools.viridian.reporting;

import static com.hartwig.hmftools.common.purple.PurpleQCStatus.FAIL_CONTAMINATION;
import static com.hartwig.hmftools.common.purple.PurpleQCStatus.FAIL_NO_TUMOR;
import static com.hartwig.hmftools.viridian.reporting.ClonalityStatus.CLONAL;
import static com.hartwig.hmftools.viridian.reporting.ClonalityStatus.NOT_CLONAL;
import static com.hartwig.hmftools.viridian.reporting.ClonalityStatus.NOT_PROFILED;
import static com.hartwig.hmftools.viridian.reporting.ClonalityStatus.PURPLE_FAIL_CONTAMINATION;
import static com.hartwig.hmftools.viridian.reporting.ClonalityStatus.PURPLE_FAIL_NO_TUMOR;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.viridian.profiling.VirusProfile;
import com.hartwig.hmftools.viridian.profiling.common.GenomeStats;
import com.hartwig.hmftools.viridian.profiling.common.SummaryStats;
import com.hartwig.hmftools.viridian.profiling.select.RepresentativeStatus;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.OncologyGroupInfo;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class VirusReporterTest
{
    private static final OncologyGroup GROUP = new OncologyGroup("HPV16");
    private static final VirusGenome GENOME = new VirusGenome("c1", 1000, "HPV type 16", GROUP);
    private static final OncologyGroupInfo REPORTABLE_INFO =
            new OncologyGroupInfo(GROUP, VirusType.HPV, VirusLikelihoodType.HIGH);
    private static final OncologyGroupInfo NOT_REPORTABLE_INFO = new OncologyGroupInfo(GROUP, null, null);

    private static final ClonalityInputs PURPLE_FIT = new ClonalityInputs(20.0, null);
    private static final ClonalityInputs NO_TUMOR = new ClonalityInputs(null, FAIL_NO_TUMOR);
    private static final ClonalityInputs CONTAMINATION = new ClonalityInputs(null, FAIL_CONTAMINATION);

    private static final boolean MEETS_INTEGRATION_MIN = true;
    private static final boolean REPORTED = true;

    // Thresholds are 1 integration and 0.5 copies per tumor cell. Copies per tumor cell is the representative's mean
    // depth over the expected depth per copy, here 20.

    @Test
    public void testReportedOnClonalityAtFloor()
    {
        VirusProfile profile = profile(10);
        assertEquals(
                reportableReport(0, !MEETS_INTEGRATION_MIN, profile, 0.5, CLONAL, REPORTED),
                VirusReporter.reportGroup(GROUP, profile, 0, REPORTABLE_INFO, PURPLE_FIT));
    }

    @Test
    public void testNotClonalBelowFloor()
    {
        VirusProfile profile = profile(9);
        assertEquals(
                reportableReport(0, !MEETS_INTEGRATION_MIN, profile, 0.45, NOT_CLONAL, !REPORTED),
                VirusReporter.reportGroup(GROUP, profile, 0, REPORTABLE_INFO, PURPLE_FIT));
    }

    @Test
    public void testPurpleFailureNamesClonalityStatus()
    {
        VirusProfile profile = profile(10);
        assertEquals(
                reportableReport(0, !MEETS_INTEGRATION_MIN, profile, null, PURPLE_FAIL_NO_TUMOR, !REPORTED),
                VirusReporter.reportGroup(GROUP, profile, 0, REPORTABLE_INFO, NO_TUMOR));
        assertEquals(
                reportableReport(0, !MEETS_INTEGRATION_MIN, profile, null, PURPLE_FAIL_CONTAMINATION, !REPORTED),
                VirusReporter.reportGroup(GROUP, profile, 0, REPORTABLE_INFO, CONTAMINATION));
    }

    @Test
    public void testReportedOnIntegrationAloneWhenNotProfiled()
    {
        assertEquals(
                reportableReport(1, MEETS_INTEGRATION_MIN, null, null, NOT_PROFILED, REPORTED),
                VirusReporter.reportGroup(GROUP, null, 1, REPORTABLE_INFO, PURPLE_FIT));
    }

    // Not being found in the reads explains the missing copies per tumor cell, whatever Purple's QC.
    @Test
    public void testNotProfiledOverridesPurpleFailure()
    {
        assertEquals(
                reportableReport(1, MEETS_INTEGRATION_MIN, null, null, NOT_PROFILED, REPORTED),
                VirusReporter.reportGroup(GROUP, null, 1, REPORTABLE_INFO, NO_TUMOR));
    }

    @Test
    public void testNotReportableNeverReportedDespiteEvidence()
    {
        VirusProfile profile = profile(10);
        VirusReport expected = new VirusReport(
                GROUP, null, null, 5, MEETS_INTEGRATION_MIN, profile, 0.5, CLONAL, !REPORTED);
        assertEquals(expected, VirusReporter.reportGroup(GROUP, profile, 5, NOT_REPORTABLE_INFO, PURPLE_FIT));
    }

    private static VirusReport reportableReport(
            int integrations, boolean meetsIntegrationMin, @Nullable VirusProfile profile,
            @Nullable Double copiesPerTumorCell, ClonalityStatus clonalityStatus, boolean isReported)
    {
        return new VirusReport(
                GROUP, VirusType.HPV, VirusLikelihoodType.HIGH, integrations, meetsIntegrationMin, profile,
                copiesPerTumorCell, clonalityStatus, isReported);
    }

    private static VirusProfile profile(int meanDepth)
    {
        SummaryStats depth = SummaryStats.from(new int[] { meanDepth });
        SummaryStats alignerScore = SummaryStats.from(new int[] { 30 });
        GenomeStats stats = new GenomeStats(GENOME, 1, 0, GENOME.length(), depth, alignerScore);
        return new VirusProfile(stats, RepresentativeStatus.PASS);
    }
}

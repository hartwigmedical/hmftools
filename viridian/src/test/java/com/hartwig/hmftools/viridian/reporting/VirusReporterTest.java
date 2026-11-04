package com.hartwig.hmftools.viridian.reporting;

import static com.hartwig.hmftools.common.purple.PurpleQCStatus.FAIL_CONTAMINATION;
import static com.hartwig.hmftools.common.purple.PurpleQCStatus.FAIL_NO_TUMOR;
import static com.hartwig.hmftools.viridian.reporting.ClonalityStatus.CLONAL;
import static com.hartwig.hmftools.viridian.reporting.ClonalityStatus.NOT_CLONAL;
import static com.hartwig.hmftools.viridian.reporting.ClonalityStatus.NOT_PROFILED;
import static com.hartwig.hmftools.viridian.reporting.ClonalityStatus.PURPLE_FAIL_CONTAMINATION;
import static com.hartwig.hmftools.viridian.reporting.ClonalityStatus.PURPLE_FAIL_NO_TUMOR;
import static com.hartwig.hmftools.viridian.reporting.ReportQcStatus.AMBIGUOUS_VIRUS;
import static com.hartwig.hmftools.viridian.reporting.ReportQcStatus.NO_ABNORMALITIES;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.viridian.profiling.VirusProfile;
import com.hartwig.hmftools.viridian.profiling.common.GenomeStats;
import com.hartwig.hmftools.viridian.profiling.common.SummaryStats;
import com.hartwig.hmftools.viridian.reference.Virus;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusGroup;
import com.hartwig.hmftools.viridian.reference.VirusGroupInfo;

import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class VirusReporterTest
{
    private static final VirusGroup GROUP = new VirusGroup("HPV");
    private static final Virus HPV16 = new Virus("HPV type 16", GROUP);
    private static final Virus HPV18 = new Virus("HPV type 18", GROUP);
    private static final VirusGenome GENOME = new VirusGenome("c1", 1000, HPV16);
    private static final VirusGenome OTHER_HPV16_GENOME = new VirusGenome("c2", 1000, HPV16);
    private static final VirusGenome HPV18_GENOME = new VirusGenome("c3", 1000, HPV18);
    private static final VirusGroupInfo REPORTABLE_INFO = new VirusGroupInfo(GROUP, VirusType.HPV);
    private static final VirusGroupInfo NOT_REPORTABLE_INFO = new VirusGroupInfo(GROUP, null);

    private static final ClonalityInputs PURPLE_FIT = new ClonalityInputs(20.0, null);
    private static final ClonalityInputs NO_TUMOR = new ClonalityInputs(null, FAIL_NO_TUMOR);
    private static final ClonalityInputs CONTAMINATION = new ClonalityInputs(null, FAIL_CONTAMINATION);

    private static final boolean REPORTED = true;

    // Thresholds are 1 integration and 0.5 copies per tumor cell. Copies per tumor cell is the representative's mean
    // depth over the expected depth per copy, here 20.

    @Test
    public void testReportedOnClonalityAtFloor()
    {
        VirusProfile profile = profile(10);
        assertEquals(
                reportableReport(HPV16, NO_ABNORMALITIES, 0, profile, 0.5, CLONAL, REPORTED),
                VirusReporter.reportGroup(profile, List.of(), REPORTABLE_INFO, PURPLE_FIT));
    }

    @Test
    public void testNotClonalBelowFloor()
    {
        VirusProfile profile = profile(9);
        assertEquals(
                reportableReport(HPV16, NO_ABNORMALITIES, 0, profile, 0.45, NOT_CLONAL, !REPORTED),
                VirusReporter.reportGroup(profile, List.of(), REPORTABLE_INFO, PURPLE_FIT));
    }

    @Test
    public void testPurpleFailureNamesClonalityStatus()
    {
        VirusProfile profile = profile(10);
        assertEquals(
                reportableReport(HPV16, NO_ABNORMALITIES, 0, profile, null, PURPLE_FAIL_NO_TUMOR, !REPORTED),
                VirusReporter.reportGroup(profile, List.of(), REPORTABLE_INFO, NO_TUMOR));
        assertEquals(
                reportableReport(HPV16, NO_ABNORMALITIES, 0, profile, null, PURPLE_FAIL_CONTAMINATION, !REPORTED),
                VirusReporter.reportGroup(profile, List.of(), REPORTABLE_INFO, CONTAMINATION));
    }

    @Test
    public void testReportedOnIntegrationAloneWhenNotProfiled()
    {
        assertEquals(
                reportableReport(HPV16, NO_ABNORMALITIES, 1, null, 0.0, NOT_PROFILED, REPORTED),
                VirusReporter.reportGroup(null, List.of(GENOME), REPORTABLE_INFO, PURPLE_FIT));
    }

    // Not being found in the reads explains the missing copies per tumor cell, whatever Purple's QC.
    @Test
    public void testNotProfiledOverridesPurpleFailure()
    {
        assertEquals(
                reportableReport(HPV16, NO_ABNORMALITIES, 1, null, null, NOT_PROFILED, REPORTED),
                VirusReporter.reportGroup(null, List.of(GENOME), REPORTABLE_INFO, NO_TUMOR));
    }

    @Test
    public void testNotReportableNeverReportedDespiteEvidence()
    {
        VirusProfile profile = profile(10);
        List<VirusGenome> integratedGenomes = List.of(GENOME, GENOME, GENOME, GENOME, GENOME);
        VirusReport expected = new VirusReport(HPV16, NO_ABNORMALITIES, null, 5, profile, 0.5, CLONAL, !REPORTED);
        assertEquals(expected, VirusReporter.reportGroup(profile, integratedGenomes, NOT_REPORTABLE_INFO, PURPLE_FIT));
    }

    // Integration inserts often land on an arbitrary genome among near-identical strains, so only the virus name counts.
    @Test
    public void testIntegrationOnOtherGenomeWithSameVirusNameIsNotAmbiguous()
    {
        VirusProfile profile = profile(10);
        assertEquals(
                reportableReport(HPV16, NO_ABNORMALITIES, 1, profile, 0.5, CLONAL, REPORTED),
                VirusReporter.reportGroup(profile, List.of(OTHER_HPV16_GENOME), REPORTABLE_INFO, PURPLE_FIT));
    }

    @Test
    public void testRepresentativeVirusNameWinsOverDisagreeingIntegrations()
    {
        VirusProfile profile = profile(10);
        assertEquals(
                reportableReport(HPV16, AMBIGUOUS_VIRUS, 2, profile, 0.5, CLONAL, REPORTED),
                VirusReporter.reportGroup(profile, List.of(HPV18_GENOME, HPV18_GENOME), REPORTABLE_INFO, PURPLE_FIT));
    }

    @Test
    public void testMostFrequentIntegrationVirusNameWhenNotProfiled()
    {
        assertEquals(
                reportableReport(HPV18, AMBIGUOUS_VIRUS, 3, null, 0.0, NOT_PROFILED, REPORTED),
                VirusReporter.reportGroup(null, List.of(HPV18_GENOME, GENOME, HPV18_GENOME), REPORTABLE_INFO, PURPLE_FIT));
    }

    @Test
    public void testTiedIntegrationVirusNamesPickAlphabeticallyFirst()
    {
        assertEquals(
                reportableReport(HPV16, AMBIGUOUS_VIRUS, 2, null, 0.0, NOT_PROFILED, REPORTED),
                VirusReporter.reportGroup(null, List.of(HPV18_GENOME, GENOME), REPORTABLE_INFO, PURPLE_FIT));
    }

    private static VirusReport reportableReport(
            Virus virus, ReportQcStatus qcStatus, int integrations, @Nullable VirusProfile profile,
            @Nullable Double copiesPerTumorCell, ClonalityStatus clonalityStatus, boolean isReported)
    {
        return new VirusReport(
                virus, qcStatus, VirusType.HPV, integrations, profile, copiesPerTumorCell, clonalityStatus,
                isReported);
    }

    private static VirusProfile profile(int meanDepth)
    {
        SummaryStats depth = SummaryStats.from(new int[] { meanDepth });
        SummaryStats alignerScore = SummaryStats.from(new int[] { 30 });
        GenomeStats stats = new GenomeStats(GENOME, 1, 0, GENOME.length(), depth, alignerScore);
        return new VirusProfile(stats, List.of());
    }
}

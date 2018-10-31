package com.hartwig.hmftools.patientreporter.copynumber;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestFactory.createTestCopyNumberBuilder;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testSequencedReportData;
import static com.hartwig.hmftools.patientreporter.copynumber.ReportingCopyNumberFilters.includeInReport;
import static com.hartwig.hmftools.patientreporter.copynumber.ReportingCopyNumberFilters.isSignificant;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;

import org.junit.Test;

public class ReportingCopyNumberFiltersTest {

    @Test
    public void canDetermineSignificantEvent() {
        assertTrue(isSignificant(1, 0));
        assertTrue(isSignificant(2, 0));
        assertFalse(isSignificant(2, 1));

        assertFalse(isSignificant(2, 2));
        assertTrue(isSignificant(2, 20));
        assertFalse(isSignificant(10, 20));
    }

    @Test
    public void canFilterForReportingCorrectly() {
        GeneModel geneModel = testSequencedReportData().panelGeneModel();
        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();

        String gene = "DOES_NOT_EXIST";
        assertFalse(geneModel.isDeletionReportable(gene));
        geneCopyNumbers.add(createTestCopyNumberBuilder().gene(gene).minCopyNumber(0).build());

        List<GeneCopyNumber> filtered = ReportingCopyNumberFilters.filterForReporting(geneCopyNumbers, geneModel);

        assertTrue(filtered.isEmpty());
    }

    @Test
    public void reportLossesCorrectly() {
        assertTrue(includeInReport(0, true, true));
        assertFalse(includeInReport(0, true, false));
        assertFalse(includeInReport(5, false, true));
    }

    @Test
    public void reportGainsCorrectly() {
        assertTrue(includeInReport(5, true, true));
        assertFalse(includeInReport(5, false, true));
        assertFalse(includeInReport(0, true, false));
    }
}

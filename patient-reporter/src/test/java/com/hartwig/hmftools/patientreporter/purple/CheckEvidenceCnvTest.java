package com.hartwig.hmftools.patientreporter.purple;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestFactory.createTestEvidenceBuilder;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.cnv.CopyNumberEvidenceAnalyzerTest;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CheckEvidenceCnvTest {

    @Test
    public void doNotFilterEvidenceOnReportableGene() {
        String gene = "GENE1";
        List<ReportableGainLoss> reportableGainsAndLosses = createTestReportableGainsAndLosses(gene);
        List<EvidenceItem> evidenceItems = createTestEvidenceItems();

        Map<ReportableGainLoss, List<EvidenceItem>> evidencePerGeneCopyNumber = Maps.newHashMap();
        evidencePerGeneCopyNumber.put(reportableGainsAndLosses.get(0), evidenceItems);

        Map<ReportableGainLoss, List<EvidenceItem>> filteredEvidenceItemMap =
                CheckEvidenceCnv.checkAndFilterForEvidenceInDriverCatalog(reportableGainsAndLosses, evidencePerGeneCopyNumber);

        assertEquals(1, filteredEvidenceItemMap.size());
    }

    @NotNull
    private static List<ReportableGainLoss> createTestReportableGainsAndLosses(@NotNull String gene) {
        return Lists.newArrayList(CopyNumberEvidenceAnalyzerTest.createTestReportableGainLossBuilder().gene(gene).build());
    }

    @NotNull
    private static List<EvidenceItem> createTestEvidenceItems() {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();

        evidenceItems.add(createTestEvidenceBuilder().build());
        evidenceItems.add(createTestEvidenceBuilder().build());
        evidenceItems.add(createTestEvidenceBuilder().build());

        return evidenceItems;
    }
}
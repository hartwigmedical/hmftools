package com.hartwig.hmftools.patientreporter.purple;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestFactory.createTestCopyNumberBuilder;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestFactory.createTestEvidenceBuilder;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestFactory.createTestReportableGainLossBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CheckEvidenceCnvTest {

    @Test
    public void canFilterOutEvidence() {
        List<ReportableGainLoss> reportableGainsAndLosses = createTestReportableGainsAndLosses("GENE1");
        List<EvidenceItem> evidenceItems = createTestEvidenceItems();

        GeneCopyNumber geneCopyNumber = createTestCopyNumberBuilder().gene("NOT_GENE1").build();

        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber = Maps.newHashMap();
        evidencePerGeneCopyNumber.put(geneCopyNumber, evidenceItems);

        Map<GeneCopyNumber, List<EvidenceItem>> filteredEvidenceItemMap =
                CheckEvidenceCnv.checkAndFilterForEvidenceInDriverCatalog(reportableGainsAndLosses, evidencePerGeneCopyNumber);

        assertTrue(filteredEvidenceItemMap.isEmpty());
    }

    @Test
    public void doNotFilterEvidenceOnReportableGene() {
        String gene = "GENE1";
        List<ReportableGainLoss> reportableGainsAndLosses = createTestReportableGainsAndLosses(gene);
        List<EvidenceItem> evidenceItems = createTestEvidenceItems();

        GeneCopyNumber geneCopyNumber = createTestCopyNumberBuilder().gene(gene).build();

        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber = Maps.newHashMap();
        evidencePerGeneCopyNumber.put(geneCopyNumber, evidenceItems);

        Map<GeneCopyNumber, List<EvidenceItem>> filteredEvidenceItemMap =
                CheckEvidenceCnv.checkAndFilterForEvidenceInDriverCatalog(reportableGainsAndLosses, evidencePerGeneCopyNumber);

        assertEquals(1, filteredEvidenceItemMap.size());
    }

    @NotNull
    private static List<ReportableGainLoss> createTestReportableGainsAndLosses(@NotNull String gene) {
        return Lists.newArrayList(createTestReportableGainLossBuilder().gene(gene).build());
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
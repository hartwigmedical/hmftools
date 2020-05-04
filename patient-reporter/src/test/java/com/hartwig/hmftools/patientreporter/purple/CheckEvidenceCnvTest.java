package com.hartwig.hmftools.patientreporter.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.EvidenceScope;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CheckEvidenceCnvTest {

    @Test
    public void canFilterOutEvidence() {
        List<ReportableGainLoss> reportableGainsAndLosses = createTestReportableGainsAndLosses("GENE1");
        List<EvidenceItem> evidenceItems = createTestEvidenceItems();

        GeneCopyNumber geneCopyNumber = createTestCopyNumberBuilder("GENE2").build();

        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber = Maps.newHashMap();
        evidencePerGeneCopyNumber.put(geneCopyNumber, evidenceItems);

        Map<GeneCopyNumber, List<EvidenceItem>> filteredEvidenceItemMap =
                CheckEvidenceCnv.checkAndFilterForEvidenceInDriverCatalog(reportableGainsAndLosses, evidencePerGeneCopyNumber);

        assertTrue(filteredEvidenceItemMap.isEmpty());
    }

    @Test
    public void doNotFilterEvidenceOnReportableGene() {
        List<ReportableGainLoss> reportableGainsAndLosses = createTestReportableGainsAndLosses("GENE1");
        List<EvidenceItem> evidenceItems = createTestEvidenceItems();

        GeneCopyNumber geneCopyNumber = createTestCopyNumberBuilder("GENE1").build();

        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber = Maps.newHashMap();
        evidencePerGeneCopyNumber.put(geneCopyNumber, evidenceItems);

        Map<GeneCopyNumber, List<EvidenceItem>> filteredEvidenceItemMap =
                CheckEvidenceCnv.checkAndFilterForEvidenceInDriverCatalog(reportableGainsAndLosses, evidencePerGeneCopyNumber);

        assertEquals(1, filteredEvidenceItemMap.size());
    }

    @NotNull
    private static List<ReportableGainLoss> createTestReportableGainsAndLosses(@NotNull String gene) {
        ReportableGainLoss gainLoss = ImmutableReportableGainLoss.builder()
                .chromosome("1")
                .chromosomeBand("band")
                .gene(gene)
                .copies(0)
                .interpretation(CopyNumberInterpretation.PARTIAL_LOSS)
                .build();

        return Lists.newArrayList(gainLoss);
    }

    @NotNull
    private static List<EvidenceItem> createTestEvidenceItems() {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();

        evidenceItems.add(createTestEvidenceBuilder().build());
        evidenceItems.add(createTestEvidenceBuilder().build());
        evidenceItems.add(createTestEvidenceBuilder().build());

        return evidenceItems;
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder createTestEvidenceBuilder() {
        return ImmutableEvidenceItem.builder()
                .event(Strings.EMPTY)
                .level(EvidenceLevel.LEVEL_A)
                .response(Strings.EMPTY)
                .reference(Strings.EMPTY)
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .drug(Strings.EMPTY)
                .drugsType(Strings.EMPTY)
                .cancerType(Strings.EMPTY)
                .isOnLabel(false);
    }

    @NotNull
    private static ImmutableGeneCopyNumber.Builder createTestCopyNumberBuilder(@NotNull String gene) {
        return ImmutableGeneCopyNumber.builder()
                .start(1)
                .end(2)
                .gene(gene)
                .chromosome("1")
                .chromosomeBand("band")
                .minRegionStart(0)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEnd(0)
                .minRegionEndSupport(SegmentSupport.NONE)
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .minRegions(1)
                .germlineHet2HomRegions(0)
                .germlineHomRegions(0)
                .somaticRegions(1)
                .minCopyNumber(0.1)
                .maxCopyNumber(0.1)
                .transcriptID("trans")
                .transcriptVersion(0)
                .minMinorAllelePloidy(0);
    }
}
package com.hartwig.hmftools.patientreporter.copynumber;

import static org.junit.Assert.*;

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

        List<ReportableGainLoss> reportableGainsAndLosses = reportableGainsAndLosses("FLT1");
        List<EvidenceItem> evidenceItems = evidenceItems();

        GeneCopyNumber geneCopyNumber = createTestCopyNumberBuilder("EGFR").build();

        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber = Maps.newHashMap();
        evidencePerGeneCopyNumber.put(geneCopyNumber, evidenceItems);

        Map<GeneCopyNumber, List<EvidenceItem>> filteredEvidenceItemMap =
                CheckEvidenceCnv.checkAndFilterForEvidenceInDriverCatalog(reportableGainsAndLosses, evidencePerGeneCopyNumber);

    }

    @Test
    public void haveNonFilterOut() {

        List<ReportableGainLoss> reportableGainsAndLosses = reportableGainsAndLosses("EGFR");
        List<EvidenceItem> evidenceItems = evidenceItems();

        GeneCopyNumber geneCopyNumber = createTestCopyNumberBuilder("EGFR").build();

        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber = Maps.newHashMap();
        evidencePerGeneCopyNumber.put(geneCopyNumber, evidenceItems);

        Map<GeneCopyNumber, List<EvidenceItem>> filteredEvidenceItemMap =
                CheckEvidenceCnv.checkAndFilterForEvidenceInDriverCatalog(reportableGainsAndLosses, evidencePerGeneCopyNumber);

    }

    @NotNull
    private static List<ReportableGainLoss> reportableGainsAndLosses(@NotNull String gene) {
        ReportableGainLoss gainLoss = ImmutableReportableGainLoss.builder()
                .chromosome("10")
                .chromosomeBand("q23.31")
                .gene(gene)
                .copies(0)
                .interpretation(CopyNumberInterpretation.PARTIAL_LOSS)
                .build();
        return Lists.newArrayList(gainLoss);
    }

    @NotNull
    private static List<EvidenceItem> evidenceItems() {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();

        ImmutableEvidenceItem.Builder offLabelBuilder = evidenceBuilder().isOnLabel(true);

        evidenceItems.add(offLabelBuilder.event("FLT1 p.Val600Glu")
                .drug("Alpelisib + Cetuximab + Encorafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:17")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.GENE_LEVEL)
                .build());

        evidenceItems.add(offLabelBuilder.event("EGFR p.Val300Glu")
                .drug("Bevacizumab")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("PTEN p.Val499Glu")
                .drug("Alpelisib + Cetuximab + Encorafenib")
                .level(EvidenceLevel.LEVEL_C)
                .response("Responsive")
                .reference("variant:17")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.GENE_LEVEL)
                .build());


        return evidenceItems;
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder evidenceBuilder() {
        return ImmutableEvidenceItem.builder().drugsType(Strings.EMPTY).cancerType(Strings.EMPTY);
    }

    @NotNull
    private static ImmutableGeneCopyNumber.Builder createTestCopyNumberBuilder(@NotNull final String gene) {
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
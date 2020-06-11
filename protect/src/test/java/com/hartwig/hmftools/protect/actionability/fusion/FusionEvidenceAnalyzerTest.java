package com.hartwig.hmftools.protect.actionability.fusion;

import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.protect.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.fusion.ImmutableReportableGeneFusion;
import com.hartwig.hmftools.common.fusion.ReportableGeneFusion;
import com.hartwig.hmftools.protect.actionability.cancertype.CancerTypeAnalyzerTestFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FusionEvidenceAnalyzerTest {

    @Test
    public void actionabilityWorksFusions() {
        ActionableFusion fusion = ImmutableActionableFusion.builder()
                .fiveGene("BCR")
                .threeGene("ABL1")
                .cancerType("Acute lymphoblastic leukemia")
                .drug("Dasatinib")
                .drugsType("BCR-ABL inhibitor")
                .level("A")
                .reference("BCR__ABL1")
                .response("Responsive")
                .source("CGI")
                .build();

        ActionablePromiscuous five = ImmutableActionablePromiscuous.builder()
                .cancerType("Acute lymphoblastic leukemia")
                .drug("FORT-1")
                .drugsType("Trial")
                .gene("FGFR1")
                .level("B")
                .reference("EXT8753 (NL62741.031.17)")
                .response("Responsive")
                .source("iclusion")
                .build();

        ActionablePromiscuous three = ImmutableActionablePromiscuous.builder()
                .cancerType("Acute lymphoblastic leukemia")
                .drug("AZD5363")
                .drugsType("AKT inhibitor")
                .gene("AKT3__.")
                .level("D")
                .reference("AKT3__.")
                .response("Responsive")
                .source("cgi")
                .build();

        FusionEvidenceAnalyzer analyzer =
                new FusionEvidenceAnalyzer((Lists.newArrayList(fusion)), Lists.newArrayList(five), Lists.newArrayList(three));

        CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzerTestFactory.buildWithOneCancerTypeMapping("Skin", "4159");

        ReportableGeneFusion reportableFusion = createTestFusionBuilder().geneStart("TMPRSS2").geneEnd("PNPLA7").build();

        assertTrue(analyzer.evidenceForFusion(reportableFusion, "Skin", cancerTypeAnalyzer).isEmpty());
    }

    @NotNull
    private static ImmutableReportableGeneFusion.Builder createTestFusionBuilder() {
        return ImmutableReportableGeneFusion.builder()
                .geneContextStart(Strings.EMPTY)
                .geneContextEnd(Strings.EMPTY)
                .geneTranscriptStart(Strings.EMPTY)
                .geneTranscriptEnd(Strings.EMPTY)
                .junctionCopyNumber(1D);
    }
}
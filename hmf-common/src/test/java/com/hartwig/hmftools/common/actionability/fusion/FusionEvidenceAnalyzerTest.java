package com.hartwig.hmftools.common.actionability.fusion;

import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzerTestFactory;
import com.hartwig.hmftools.common.variant.structural.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.variant.structural.linx.FusionPhasedType;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;

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

        LinxFusion reportableFusion = createTestFusionBuilder().geneStart("TMPRSS2").geneEnd("PNPLA7").build();

        assertTrue(analyzer.evidenceForFusion(reportableFusion, "Skin", cancerTypeAnalyzer).isEmpty());
    }

    @NotNull
    private static ImmutableLinxFusion.Builder createTestFusionBuilder() {
        return ImmutableLinxFusion.builder()
                .fivePrimeBreakendId(1)
                .threePrimeBreakendId(2)
                .name(Strings.EMPTY)
                .reported(true)
                .reportedType(Strings.EMPTY)
                .phased(FusionPhasedType.SKIPPED_EXONS)
                .likelihood(FusionLikelihoodType.HIGH)
                .chainLength(1)
                .chainLinks(1)
                .chainTerminated(true)
                .domainsKept(Strings.EMPTY)
                .domainsLost(Strings.EMPTY)
                .skippedExonsUp(2)
                .skippedExonsDown(4)
                .fusedExonUp(6)
                .fusedExonDown(7)
                .geneStart(Strings.EMPTY)
                .geneContextStart(Strings.EMPTY)
                .geneTranscriptStart(Strings.EMPTY)
                .geneEnd(Strings.EMPTY)
                .geneContextEnd(Strings.EMPTY)
                .geneTranscriptEnd(Strings.EMPTY)
                .junctionCopyNumber(1D);
    }
}
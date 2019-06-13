package com.hartwig.hmftools.common.actionability.fusion;

import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzerTestFactory;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableSimpleGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.SimpleGeneFusion;

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

        SimpleGeneFusion simpleGeneFusion = ImmutableSimpleGeneFusion.builder()
                .fiveGene("TMPRSS2")
                .threeGene("PNPLA7")
                .build();

        assertTrue(analyzer.evidenceForFusion(simpleGeneFusion, "Skin", cancerTypeAnalyzer).isEmpty());
    }
}
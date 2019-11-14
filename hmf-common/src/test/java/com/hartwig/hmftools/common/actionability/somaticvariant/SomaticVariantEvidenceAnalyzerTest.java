package com.hartwig.hmftools.common.actionability.somaticvariant;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzerTestFactory;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariant;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class SomaticVariantEvidenceAnalyzerTest {

    @Test
    public void actionabilityWorksVariants() {
        ActionableSomaticVariant actionableSomaticVariant = ImmutableActionableSomaticVariant.builder()
                .gene("BRAF")
                .chromosome("X")
                .position(1234)
                .ref("C")
                .alt("T")
                .source("civic")
                .reference("BRAF 600E")
                .drug("Dabrafenib")
                .drugsType("BRAF inhibitor")
                .cancerType("Skin Melanoma")
                .level("A")
                .response("Responsive")
                .build();

        ActionableRange actionableRange = ImmutableActionableRange.builder()
                .gene("BRAF")
                .chromosome("7")
                .start(10)
                .end(1500)
                .source("oncoKB")
                .reference("BRAF Oncogenic Mutations")
                .drug("Cetuximab")
                .drugsType("EGFR mAb inhibitor")
                .cancerType("Skin Melanoma")
                .level("A")
                .response("Resistant")
                .build();

        SomaticVariantEvidenceAnalyzer analyzer =
                new SomaticVariantEvidenceAnalyzer(Lists.newArrayList(actionableSomaticVariant), Lists.newArrayList(actionableRange));

        CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzerTestFactory.buildWithOneCancerTypeMapping("Skin Melanoma", "4159");

        ReportableVariant variantWithNoCodingEffect = ImmutableReportableVariant.builder()
                .gene("BRAF")
                .chromosome("7")
                .position(100)
                .ref("C")
                .alt("T")
                .totalReadCount(0)
                .alleleReadCount(0)
                .gDNA(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.NONE)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .hotspot(Hotspot.NON_HOTSPOT)
                .biallelic(false)
                .clonalLikelihood(1)
                .totalPloidy(0D)
                .allelePloidy(0D)
                .notifyClinicalGeneticist(false)
                .build();

        assertEquals(0, analyzer.evidenceForVariant(variantWithNoCodingEffect, "Skin", cancerTypeAnalyzer).size());

        ReportableVariant variantWithCodingEffect = ImmutableReportableVariant.builder()
                .gene("BRAF")
                .chromosome("7")
                .position(100)
                .ref("C")
                .alt("T")
                .totalReadCount(0)
                .alleleReadCount(0)
                .gDNA(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .hotspot(Hotspot.NON_HOTSPOT)
                .biallelic(false)
                .clonalLikelihood(1)
                .totalPloidy(0D)
                .allelePloidy(0D)
                .notifyClinicalGeneticist(false)
                .build();

        assertEquals(1, analyzer.evidenceForVariant(variantWithCodingEffect, "Skin", cancerTypeAnalyzer).size());
    }
}
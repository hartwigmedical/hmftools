package com.hartwig.hmftools.common.actionability.variant;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzerTestFactory;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class VariantEvidenceAnalyzerTest {

    @Test
    public void actionabilityWorksVariants() {
        ActionableVariant actionableVariant = ImmutableActionableVariant.builder()
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

        VariantEvidenceAnalyzer analyzer =
                new VariantEvidenceAnalyzer(Lists.newArrayList(actionableVariant), Lists.newArrayList(actionableRange));

        CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzerTestFactory.buildWithOneCancerTypeMapping("Skin Melanoma", "4159");

        SomaticVariant variantWithNoCodingEffect = SomaticVariantTestBuilderFactory.create()
                .gene("BRAF")
                .chromosome("7")
                .position(100)
                .canonicalCodingEffect(CodingEffect.NONE)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .build();

        assertEquals(0, analyzer.evidenceForVariant(variantWithNoCodingEffect, "Skin", cancerTypeAnalyzer).size());

        SomaticVariant variantWithCodingEffect = SomaticVariantTestBuilderFactory.create()
                .gene("BRAF")
                .chromosome("7")
                .position(100)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .build();

        assertEquals(1, analyzer.evidenceForVariant(variantWithCodingEffect, "Skin", cancerTypeAnalyzer).size());
    }
}
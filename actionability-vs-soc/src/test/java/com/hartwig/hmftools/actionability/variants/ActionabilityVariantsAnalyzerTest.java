package com.hartwig.hmftools.actionability.variants;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.actionability.cancerTypeMapping.CancerTypeAnalyzer;
import com.hartwig.hmftools.actionability.cancerTypeMapping.CancerTypeReading;
import com.hartwig.hmftools.actionability.cancerTypeMapping.ImmutableCancerTypeReading;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ActionabilityVariantsAnalyzerTest {

    @Test
    public void actionabilityWorksVariants() {
        ActionabilityVariant actionabilityVariant = ImmutableActionabilityVariant.builder()
                .gene("BRAF")
                .chromosome("7")
                .position(100)
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

        ActionabilityRange actionabilityRange = ImmutableActionabilityRange.builder()
                .gene("BRAF")
                .chromosome("7")
                .start(10)
                .end(1500)
                .source("oncoKB")
                .reference("NRAS Oncogenic Mutations")
                .drug("Cetuximab")
                .drugsType("EGFR mAb inhibitor")
                .cancerType("Skin Melanoma")
                .level("A")
                .response("Resistant")
                .build();

        CancerTypeReading reading = ImmutableCancerTypeReading.builder().doidSet("4159").cancerType("Skin").build();

        ActionabilityVariantsAnalyzer analyzer =
                new ActionabilityVariantsAnalyzer(Lists.newArrayList(actionabilityVariant), Lists.newArrayList(actionabilityRange));

        CancerTypeAnalyzer cancerType = new CancerTypeAnalyzer(Lists.newArrayList(reading));

        SomaticVariant variant = ImmutableSomaticVariantImpl.builder()
                .chromosome("7")
                .position(100)
                .ref("C")
                .alt("T")
                .type(VariantType.UNDEFINED)
                .filter(Strings.EMPTY)
                .gene("BRAF")
                .genesEffected(0)
                .worstEffect(Strings.EMPTY)
                .worstEffectTranscript(Strings.EMPTY)
                .worstCodingEffect(CodingEffect.NONE)
                .totalReadCount(0)
                .alleleReadCount(0)
                .canonicalEffect(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .hotspot(Hotspot.NON_HOTSPOT)
                .mappability(0D)
                .build();

        assertEquals(true, analyzer.actionableVariants(variant, cancerType, "4159", "Skin"));
        assertEquals(false, analyzer.actionableVariants(variant, cancerType, "4159", "Breast"));
        assertEquals(true, analyzer.actionableRange(variant, cancerType, "4159", "Skin"));
        assertEquals(false, analyzer.actionableRange(variant, cancerType, "4159", "Kidney"));
    }
}
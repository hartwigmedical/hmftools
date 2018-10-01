package com.hartwig.hmftools.actionability.variants;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

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
        ActionabilityVariantsSOC variantsSOC = ImmutableActionabilityVariantsSOC.builder()
                .gene("BRAF")
                .chromosome("7")
                .position("100")
                .ref("C")
                .alt("T")
                .source("civic")
                .reference("BRAF 600E")
                .drug("Dabrafenib")
                .drugsType("BRAF inhibitor")
                .cancerType("Skin Melanoma")
                .levelSource("1")
                .levelHmf("A")
                .evidenceType("Predictive")
                .significanceSource("Responsive")
                .hmfResponse("Responsive")
                .build();

        ActionabilityRanges variantsRanges = ImmutableActionabilityRanges.builder()
                .gene("BRAF")
                .mutationTranscript("ENST00000357654")
                .chromosome("7")
                .start(Integer.toString(10))
                .stop(Integer.toString(1500))
                .geneTranscript("ENST00000256078")
                .source("oncoKB")
                .reference("NRAS Oncogenic Mutations")
                .drugsName("Cetuximab")
                .drugsType("EGFR mAb inhibitor")
                .cancerType("Skin Melanoma")
                .levelSource("1")
                .hmfLevel("A")
                .evidenceType("Predictive")
                .significanceSource("Resistant")
                .hmfResponse("Resistant")
                .build();

        CancerTypeReading reading = ImmutableCancerTypeReading.builder()
                .doidSet("4159")
                .cancerType("Skin")
                .build();

        ActionabilityVariantsAnalyzer var = new ActionabilityVariantsAnalyzer(Lists.newArrayList(variantsSOC), Lists.newArrayList(variantsRanges));

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


        assertTrue(var.actionableVariants(variant, cancerType, "4159", "Skin"));
      //  assertFalse(var.actionableVariants(variant, cancerType, "4159", "Stomach"));

        assertTrue(var.actionableRange(variant, cancerType, "4159", "Skin"));
    //   assertFalse(var.actionableRange(variant, cancerType, "4159", "Stomach"));
    }
}
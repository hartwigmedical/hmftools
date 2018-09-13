package com.hartwig.hmftools.actionability.variants;

import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ActionabilityVariantsAnalyzerTest {

    @Test
    public void actionabilityWorks() {

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

        ActionabilityVariantsAnalyzer var = new ActionabilityVariantsAnalyzer(Lists.newArrayList(variantsSOC));

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
                .hotspot(false)
                .mappability(0D)
                .build();

        assertTrue(var.actionable(variant, "Skin Melanoma", 1));
    }
}
package com.hartwig.hmftools.protect.actionability.variant;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.protect.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.protect.actionability.cancertype.CancerTypeAnalyzerTestFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
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

        SomaticVariant variantWithNoCodingEffect = createSomaticvariant()
                .gene("BRAF")
                .chromosome("7")
                .position(100)
                .canonicalCodingEffect(CodingEffect.NONE)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .build();

        assertEquals(0, analyzer.evidenceForVariant(variantWithNoCodingEffect, "Skin", cancerTypeAnalyzer).size());

        SomaticVariant variantWithCodingEffect = createSomaticvariant()
                .gene("BRAF")
                .chromosome("7")
                .position(100)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .build();

        assertEquals(1, analyzer.evidenceForVariant(variantWithCodingEffect, "Skin", cancerTypeAnalyzer).size());
    }

    @NotNull
    public static ImmutableSomaticVariantImpl.Builder createSomaticvariant() {
        return ImmutableSomaticVariantImpl.builder()
                .qual(100)
                .chromosome(Strings.EMPTY)
                .position(0L)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .type(VariantType.UNDEFINED)
                .filter(Strings.EMPTY)
                .totalReadCount(0)
                .alleleReadCount(0)
                .gene(Strings.EMPTY)
                .genesAffected(0)
                .worstEffect(Strings.EMPTY)
                .worstCodingEffect(CodingEffect.NONE)
                .worstEffectTranscript(Strings.EMPTY)
                .canonicalEffect(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .hotspot(Hotspot.NON_HOTSPOT)
                .recovered(false)
                .adjustedCopyNumber(0D)
                .adjustedVAF(0D)
                .minorAlleleCopyNumber(0D)
                .germlineStatus(GermlineStatus.UNKNOWN)
                .variantCopyNumber(0)
                .biallelic(false)
                .kataegis(Strings.EMPTY)
                .trinucleotideContext(Strings.EMPTY)
                .highConfidenceRegion(false)
                .microhomology(Strings.EMPTY)
                .repeatSequence(Strings.EMPTY)
                .repeatCount(0)
                .subclonalLikelihood(0)
                .tier(VariantTier.UNKNOWN)
                .mappability(0D);
    }
}
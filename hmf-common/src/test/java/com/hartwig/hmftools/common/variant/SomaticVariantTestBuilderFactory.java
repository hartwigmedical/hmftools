package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.purple.region.GermlineStatus;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class SomaticVariantTestBuilderFactory {

    private SomaticVariantTestBuilderFactory() {
    }

    @NotNull
    public static ImmutableSomaticVariantImpl.Builder create() {
        return ImmutableSomaticVariantImpl.builder()
                .chromosome(Strings.EMPTY)
                .position(0L)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .type(VariantType.UNDEFINED)
                .filter(Strings.EMPTY)
                .totalReadCount(0)
                .alleleReadCount(0)
                .gene(Strings.EMPTY)
                .genesEffected(0)
                .worstEffect(Strings.EMPTY)
                .worstCodingEffect(CodingEffect.NONE)
                .worstEffectTranscript(Strings.EMPTY)
                .canonicalEffect(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .hotspot(false)
                .mappability(0D);
    }

    @NotNull
    public static ImmutableEnrichedSomaticVariant.Builder createEnriched() {
        return ImmutableEnrichedSomaticVariant.builder().from(create().build())
                .trinucleotideContext(Strings.EMPTY)
                .highConfidenceRegion(false)
                .microhomology(Strings.EMPTY)
                .repeatSequence(Strings.EMPTY)
                .repeatCount(0)
                .clonality(Clonality.UNKNOWN)
                .adjustedCopyNumber(0D)
                .adjustedVAF(0D)
                .ploidy(0D)
                .minorAllelePloidy(0D)
                .biallelic(false)
                .germlineStatus(GermlineStatus.UNKNOWN);
    }
}





















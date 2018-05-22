package com.hartwig.hmftools.common.variant;

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
                .gene(Strings.EMPTY)
                .genesEffected(0)
                .worstEffect(Strings.EMPTY)
                .worstEffectTranscript(Strings.EMPTY)
                .worstCodingEffect(CodingEffect.NONE)
                .totalReadCount(0)
                .alleleReadCount(0)
                .hotspot(false)
                .mappability(0D);
    }
}

package com.hartwig.hmftools.common.protect.variant;

import com.hartwig.hmftools.common.variant.CodingEffect;

import org.jetbrains.annotations.NotNull;

public final class OtherEffectsInterpreter {

    private static final String DELIMITER = "\\|";

    private static final int TRANSCRIPT_INDEX = 0;
    private static final int HGVS_CODING_IMPACT_INDEX = 1;
    private static final int HGVS_PROTEIN_IMPACT_INDEX = 2;
    private static final int EFFECT_INDEX = 3;
    private static final int CODING_EFFECT_INDEX = 4;

    private OtherEffectsInterpreter() {
    }

    @NotNull
    public static String transcript(@NotNull String otherReportedEffects) {
        return otherReportedEffects.split(DELIMITER)[TRANSCRIPT_INDEX];
    }

    @NotNull
    public static String hgvsCodingImpact(@NotNull String otherReportedEffects) {
        return otherReportedEffects.split(DELIMITER)[HGVS_CODING_IMPACT_INDEX];
    }

    @NotNull
    public static String hgvsProteinImpact(@NotNull String otherReportedEffects) {
        return otherReportedEffects.split(DELIMITER)[HGVS_PROTEIN_IMPACT_INDEX];
    }

    @NotNull
    public static String effect(@NotNull String otherReportedEffects) {
        return otherReportedEffects.split(DELIMITER)[EFFECT_INDEX];
    }

    @NotNull
    public static CodingEffect codingEffect(@NotNull String otherReportedEffects) {
        return CodingEffect.valueOf(otherReportedEffects.split(DELIMITER)[CODING_EFFECT_INDEX]);
    }
}

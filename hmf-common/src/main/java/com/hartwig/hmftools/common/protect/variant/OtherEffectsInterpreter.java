package com.hartwig.hmftools.common.protect.variant;

import static com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo.VAR_IMPACT_OTHER_REPORT_DELIM;
import static com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo.VAR_IMPACT_OTHER_REPORT_ITEM_DELIM;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class OtherEffectsInterpreter {

    private static final String IMPACT_SPLITTER = VAR_IMPACT_OTHER_REPORT_DELIM;
    private static final String FIELD_SPLITTER = "\\" + VAR_IMPACT_OTHER_REPORT_ITEM_DELIM;

    private static final int TRANSCRIPT_INDEX = 0;
    private static final int HGVS_CODING_IMPACT_INDEX = 1;
    private static final int HGVS_PROTEIN_IMPACT_INDEX = 2;
    private static final int EFFECT_INDEX = 3;
    private static final int CODING_EFFECT_INDEX = 4;

    private OtherEffectsInterpreter() {
    }

    @NotNull
    public static String transcript(@NotNull String otherReportedEffects) {
        if (otherReportedEffects.isEmpty()) {
            return Strings.EMPTY;
        }
        return first(otherReportedEffects).split(FIELD_SPLITTER)[TRANSCRIPT_INDEX];
    }

    @NotNull
    public static String hgvsCodingImpact(@NotNull String otherReportedEffects) {
        if (otherReportedEffects.isEmpty()) {
            return Strings.EMPTY;
        }
        return first(otherReportedEffects).split(FIELD_SPLITTER)[HGVS_CODING_IMPACT_INDEX];
    }

    @NotNull
    public static String hgvsProteinImpact(@NotNull String otherReportedEffects) {
        if (otherReportedEffects.isEmpty()) {
            return Strings.EMPTY;
        }
        return first(otherReportedEffects).split(FIELD_SPLITTER)[HGVS_PROTEIN_IMPACT_INDEX];
    }

    @NotNull
    public static String effect(@NotNull String otherReportedEffects) {
        if (otherReportedEffects.isEmpty()) {
            return Strings.EMPTY;
        }
        return first(otherReportedEffects).split(FIELD_SPLITTER)[EFFECT_INDEX];
    }

    @NotNull
    public static CodingEffect codingEffect(@NotNull String otherReportedEffects) {
        if (otherReportedEffects.isEmpty()) {
            return CodingEffect.UNDEFINED;
        }
        return CodingEffect.valueOf(first(otherReportedEffects).split(FIELD_SPLITTER)[CODING_EFFECT_INDEX]);
    }

    @NotNull
    private static String first(@NotNull String otherReportedEffects) {
        // This is a little ugly and not full-proof but if there is one entry that contains IMPACT_SPLITTER it will cause issues that are
        // prevented in case there is just one entry.
        if (otherReportedEffects.split(FIELD_SPLITTER).length > 5) {
            return otherReportedEffects.split(IMPACT_SPLITTER)[0];
        }

        return otherReportedEffects;
    }
}

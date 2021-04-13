package com.hartwig.hmftools.ckb.classification;

import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

import org.jetbrains.annotations.NotNull;

public final class EventAndGeneExtractor {

    private EventAndGeneExtractor() {
    }

    @NotNull
    public static String extractGene(@NotNull Variant variant) {
        if (isFusion(variant) && !isPromiscuousFusion(variant)) {
            return variant.fullName();
        } else {
            return variant.gene().geneSymbol();
        }
    }

    @NotNull
    public static String extractEvent(@NotNull Variant variant) {
        if (isPromiscuousFusion(variant)) {
            return "fusion promiscuous";
        } else if (isFusion(variant)) {
            return removeAllSpaces(variant.variant());
        } else if (variant.variant().contains("exon") && !variant.variant().contains("exon ")) {
            // Some exon variants do not contain a space after the "exon" keyword
            return variant.variant().replace("exon", "exon ");
        } else {
            return variant.variant();
        }
    }

    private static boolean isFusion(@NotNull Variant variant) {
        return variant.impact() != null && variant.impact().equals("fusion");
    }

    private static boolean isPromiscuousFusion(@NotNull Variant variant) {
        return isFusion(variant) && variant.variant().equals("fusion");
    }

    @NotNull
    private static String removeAllSpaces(@NotNull String value) {
        return value.replaceAll("\\s+", "");
    }
}

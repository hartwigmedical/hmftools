package com.hartwig.hmftools.serve.extraction.util;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class KeyFormatter {

    private KeyFormatter() {
    }

    @NotNull
    public static String toProteinKey(@NotNull String gene, @Nullable String transcript, @NotNull String proteinAnnotation) {
        String formattedProteinAnnotation = !proteinAnnotation.isEmpty() ? "p." + proteinAnnotation : "-";
        return gene + "|" + transcript + "|" + formattedProteinAnnotation;
    }

    @NotNull
    public static String toExonKey(@NotNull String gene, @Nullable String transcript, int exonIndex) {
        return gene + "|" + transcript + "|" + exonIndex;
    }

    @NotNull
    public static String toCodonKey(@NotNull String gene, @Nullable String transcript, int codonIndex) {
        return gene + "|" + transcript + "|" + codonIndex;
    }
}

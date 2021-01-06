package com.hartwig.hmftools.serve.extraction.hotspot;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ProteinKeyFormatter {

    private ProteinKeyFormatter() {
    }

    @NotNull
    public static String toProteinKey(@NotNull String gene, @Nullable String transcript, @NotNull String proteinAnnotation) {
        String formattedProteinAnnotation = !proteinAnnotation.isEmpty() ? "p." + proteinAnnotation : "-";
        return gene + "|" + transcript + "|" + formattedProteinAnnotation;
    }

    @NotNull
    public static String toExonKey(@NotNull String gene, @Nullable String transcript, @NotNull String exonIndex) {
        String formattedExonIndex = !exonIndex.isEmpty() ? exonIndex : "-";
        return gene + "|" + transcript + "|" + formattedExonIndex;
    }
}

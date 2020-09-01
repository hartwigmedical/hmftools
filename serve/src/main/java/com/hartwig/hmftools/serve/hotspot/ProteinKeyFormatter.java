package com.hartwig.hmftools.serve.hotspot;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ProteinKeyFormatter {

    private ProteinKeyFormatter() {
    }

    @NotNull
    public static String toProteinKey(@NotNull HotspotAnnotation annotation) {
        return toProteinKey(annotation.gene(), annotation.transcript(), annotation.proteinAnnotation());
    }

    @NotNull
    public static String toProteinKey(@NotNull String gene, @Nullable String transcript, @NotNull String proteinAnnotation) {
        String formattedProteinAnnotation = !proteinAnnotation.isEmpty() ? "p." + proteinAnnotation : "-";
        return gene + "|" + transcript + "|" + formattedProteinAnnotation;
    }
}

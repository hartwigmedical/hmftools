package com.hartwig.hmftools.serve.util;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ProteinKeyFormatter {

    private ProteinKeyFormatter() {
    }

    @NotNull
    public static String toProteinKey(@NotNull String gene, @Nullable String transcript, @NotNull String proteinAnnotation) {
        return gene + "|" + transcript + "|p." + proteinAnnotation;
    }
}

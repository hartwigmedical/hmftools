package com.hartwig.hmftools.knowledgebasegenerator.transvar;

import org.jetbrains.annotations.NotNull;

final class TransvarConverter {

    private static final String FIELD_DELIMITER = "\t";

    private TransvarConverter() {
    }

    @NotNull
    static TransvarRecord toTransvarRecord(@NotNull String transvarLine) {
        String[] fields = transvarLine.split(FIELD_DELIMITER);

        return null;

    }
}

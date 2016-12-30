package com.hartwig.healthchecker.runners;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import org.jetbrains.annotations.NotNull;

public enum CheckType {
    METADATA,
    MAPPING,
    PRESTATS,
    KINSHIP,
    INSERT_SIZE,
    SUMMARY_METRICS,
    WGS_METRICS,
    SOMATIC,
    SLICED,
    REALIGNER,
    GERMLINE,
    COPYNUMBER;

    @NotNull
    public static Optional<CheckType> getByCategory(@NotNull final String typeToCheck) {
        final List<CheckType> types = Arrays.asList(CheckType.values());

        return types.stream().filter(type -> type.toString().equalsIgnoreCase(typeToCheck)).findFirst();
    }
}

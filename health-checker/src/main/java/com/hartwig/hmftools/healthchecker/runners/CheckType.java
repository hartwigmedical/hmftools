package com.hartwig.hmftools.healthchecker.runners;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import org.jetbrains.annotations.NotNull;

public enum CheckType {
    KINSHIP,
    COVERAGE,
    SOMATIC_VARIANTS;

    @NotNull
    public static Optional<CheckType> getByCategory(@NotNull final String typeToCheck) {
        final List<CheckType> types = Arrays.asList(CheckType.values());

        return types.stream().filter(type -> type.toString().equalsIgnoreCase(typeToCheck)).findFirst();
    }
}

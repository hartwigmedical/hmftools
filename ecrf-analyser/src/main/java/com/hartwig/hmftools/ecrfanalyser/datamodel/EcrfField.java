package com.hartwig.hmftools.ecrfanalyser.datamodel;

import java.util.Map;

import org.jetbrains.annotations.NotNull;

public class EcrfField {
    @NotNull
    private final String category;
    @NotNull
    private final String fieldName;
    @NotNull
    private final String description;
    @NotNull
    private final Map<Integer, String> valuesMap;

    public EcrfField(@NotNull final String category, @NotNull final String fieldName,
            @NotNull final String description, @NotNull final Map<Integer, String> valuesMap) {
        this.category = category;
        this.fieldName = fieldName;
        this.description = description;
        this.valuesMap = valuesMap;
    }
}

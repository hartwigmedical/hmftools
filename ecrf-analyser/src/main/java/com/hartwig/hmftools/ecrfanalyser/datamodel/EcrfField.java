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
    private final Map<Integer, String> values;

    public EcrfField(@NotNull final String category, @NotNull final String fieldName,
            @NotNull final String description, @NotNull final Map<Integer, String> values) {
        this.category = category;
        this.fieldName = fieldName;
        this.description = description;
        this.values = values;
    }

    @NotNull
    public String category() {
        return category;
    }

    @NotNull
    public String fieldName() {
        return fieldName;
    }

    @NotNull
    public String description() {
        return description;
    }

    @NotNull
    public Map<Integer, String> values() {
        return values;
    }

    public boolean isFreeText() {
        return values.size() == 0;
    }
}

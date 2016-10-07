package com.hartwig.hmftools.ecrfanalyser.datamodel;

import java.util.Map;

import org.jetbrains.annotations.NotNull;

public class EcrfField implements Comparable<EcrfField> {

    private static final String BIRTH_DATE_IDENTIFIER = "BIRTHDTC";

    @NotNull
    private final String name;
    @NotNull
    private final String description;
    @NotNull
    private final Map<Integer, String> values;

    public EcrfField(@NotNull final String name, @NotNull final String description,
            @NotNull final Map<Integer, String> values) {
        this.name = name;
        this.description = description;
        this.values = values;
    }

    @NotNull
    public String name() {
        return name;
    }

    @NotNull
    public String description() {
        return description;
    }

    @NotNull
    public Map<Integer, String> values() {
        return values;
    }

    @NotNull
    public String resolveValue(@NotNull String ecrfValue) throws EcrfResolveException {
        String value;
        if (values.size() > 0 && ecrfValue.length() > 0) {
            if (isInteger(ecrfValue)) {
                value = values.get(Integer.valueOf(ecrfValue));
                if (value == null) {
                    throw new EcrfResolveException("Could not find value in dropdown for " + name + ": " + ecrfValue);
                }
            } else {
                throw new EcrfResolveException(
                        "Could not convert value from list to integer for " + name + ": " + ecrfValue);
            }
        } else {
            value = ecrfValue;
        }

        if (name.contains(BIRTH_DATE_IDENTIFIER) && value.length() > 0) {
            if (value.length() < 4) {
                throw new EcrfResolveException("Could not convert " + name() + ": " + value);
            } else {
                value = value.substring(0, 4) + "-01-01";
            }
        }

        return value;
    }

    private static boolean isInteger(@NotNull String integerString) {
        try {
            return Integer.valueOf(integerString) != null;
        } catch (NumberFormatException exception) {
            return false;
        }
    }

    public int compareTo(@NotNull EcrfField other) {
        return name.compareTo(other.name);
    }
}

package com.hartwig.hmftools.common.purple;

import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public enum PurpleQCStatus {
    PASS,

    WARN_DELETED_GENES,
    WARN_HIGH_COPY_NUMBER_NOISE,
    WARN_GENDER_MISMATCH,
    WARN_LOW_PURITY,

    FAIL_CONTAMINATION,
    FAIL_NO_TUMOR;

    @NotNull
    public static String toString(Set<PurpleQCStatus> status) {
        return status.stream().map(Enum::toString).collect(Collectors.joining(","));
    }

    @NotNull
    public static Set<PurpleQCStatus> fromString(String line) {
        return Arrays.stream(line.split(",")).map(PurpleQCStatus::fromSingle).collect(Collectors.toSet());
    }

    @NotNull
    private static PurpleQCStatus fromSingle(String value) {
        switch (value) {
            case "FAIL_SEGMENT": return WARN_HIGH_COPY_NUMBER_NOISE;
            case "FAIL_GENDER": return WARN_GENDER_MISMATCH;
            case "FAIL_DELETED_GENES": return WARN_DELETED_GENES;
        }

        return PurpleQCStatus.valueOf(value);
    }

}
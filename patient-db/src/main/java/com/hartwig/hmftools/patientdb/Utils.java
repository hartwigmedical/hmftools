package com.hartwig.hmftools.patientdb;

import java.time.LocalDate;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class Utils {

    private Utils() {
    }

    @Nullable
    public static java.sql.Date toSQLDate(@Nullable final LocalDate date) {
        return date != null ? java.sql.Date.valueOf(date) : null;
    }

    static boolean anyNull(@NotNull final Object... arguments) {
        for (final Object object : arguments) {
            if (object == null) {
                return true;
            }
        }
        return false;
    }

    @NotNull
    static String getPatientId(@NotNull final String runName) {
        final String[] names = runName.split("_");
        return names[4];
    }
}

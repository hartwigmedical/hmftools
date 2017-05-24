package com.hartwig.hmftools.patientdb;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class Utils {

    private Utils() {
    }

    @Nullable
    static java.sql.Date toSQLDate(@Nullable final LocalDate date) {
        return date != null ? java.sql.Date.valueOf(date) : null;
    }

    public static boolean anyNotNull(@NotNull final Object... arguments) {
        for (final Object object : arguments) {
            if (object != null) {
                return true;
            }
        }
        return false;
    }

    static boolean anyNull(@NotNull final Object... arguments) {
        for (final Object object : arguments) {
            if (object == null) {
                return true;
            }
        }
        return false;
    }
}

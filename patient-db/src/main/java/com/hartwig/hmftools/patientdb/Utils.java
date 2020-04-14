package com.hartwig.hmftools.patientdb;

import java.time.LocalDate;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class Utils {

    private static final Logger LOGGER = LogManager.getLogger(Utils.class);

    private Utils() {
    }

    @Nullable
    public static java.sql.Date toSQLDate(@Nullable LocalDate date) {
        return date != null ? java.sql.Date.valueOf(date) : null;
    }

    static boolean anyNull(@NotNull Object... arguments) {
        for (Object object : arguments) {
            if (object == null) {
                return true;
            }
        }
        return false;
    }

    @NotNull
    public static String capitalize(@NotNull String string) {
        if (string.isEmpty()) {
            return string;
        } else {
            return string.toUpperCase().substring(0, 1) + string.substring(1);
        }
    }

    @NotNull
    static String extractPatientIdentifier(@NotNull String setName) {
        String[] names = setName.split("_");
        if (names.length < 5) {
            LOGGER.error("Run name {} had less than 5 parts after splitting on _", setName);
            return Strings.EMPTY;
        }
        return names[4];
    }
}

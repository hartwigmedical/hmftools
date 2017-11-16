package com.hartwig.hmftools.patientdb;

import java.time.LocalDate;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.context.RunContext;

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
    private static String getPatientId(@NotNull final String runName) {
        final String[] names = runName.split("_");
        if (names.length < 5) {
            LOGGER.error("run name {} had less than 5 parts after splitting on _", runName);
            return Strings.EMPTY;
        }
        return names[4];
    }

    @NotNull
    public static String capitalize(@NotNull final String string) {
        if (string.isEmpty()) {
            return string;
        } else {
            return string.toUpperCase().substring(0, 1) + string.substring(1);
        }
    }

    @NotNull
    static Set<String> sequencedPatientIds(@NotNull final List<RunContext> runContexts) {
        return runContexts.stream()
                .map(runContext -> getPatientId(runContext.setName()))
                .filter(patientId -> !patientId.isEmpty())
                .collect(Collectors.toSet());
    }
}

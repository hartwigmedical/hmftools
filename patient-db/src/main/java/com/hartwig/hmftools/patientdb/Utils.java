package com.hartwig.hmftools.patientdb;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class Utils {
    private static final Logger LOGGER = LogManager.getLogger(PatientDbRunner.class);

    static int getMaxLength(@NotNull Collection<List<String>> lists, @NotNull String warnMessage) {
        final Iterator<List<String>> listIterator = lists.iterator();
        int maxSize = -1;
        while (listIterator.hasNext()) {
            List<String> currentList = listIterator.next();
            if (currentList != null) {
                if (maxSize < currentList.size()) {
                    if (maxSize != -1) {
                        LOGGER.warn(warnMessage);
                    }
                    maxSize = currentList.size();
                }
            }
        }
        return maxSize;
    }

    @Nullable
    static LocalDate getDate(@Nullable String dateFieldValue, @NotNull DateTimeFormatter dateFormatter) {
        try {
            return LocalDate.parse(dateFieldValue, dateFormatter);
        } catch (NullPointerException | DateTimeParseException e) {
            return null;
        }
    }

    @Nullable
    static String getElemAtIndex(@Nullable List<String> list, int index) {
        try {
            return list.get(index);
        } catch (NullPointerException | IndexOutOfBoundsException e) {
            return null;
        }
    }

    @Nullable
    static java.sql.Date toSQLDate(@Nullable LocalDate date) {
        try {
            return java.sql.Date.valueOf(date);
        } catch (NullPointerException e) {
            return null;
        }
    }
}

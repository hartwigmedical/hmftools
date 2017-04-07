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

public class Utils {
    private static final Logger LOGGER = LogManager.getLogger(Utils.class);

    public static int getMaxLength(@NotNull final Collection<List<?>> lists, @NotNull final String warnMessage) {
        final Iterator<List<?>> listIterator = lists.iterator();
        int maxSize = -1;
        while (listIterator.hasNext()) {
            List<?> currentList = listIterator.next();
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
    public static LocalDate getDate(@Nullable final String dateFieldValue,
            @NotNull final DateTimeFormatter dateFormatter) {
        if (dateFieldValue == null) {
            return null;
        }

        try {
            return LocalDate.parse(dateFieldValue, dateFormatter);
        } catch (DateTimeParseException e) {
            return null;
        }
    }

    @Nullable
    public static <T> T getElemAtIndex(@Nullable final List<T> list, final int index) {
        if (list == null) {
            return null;
        }

        try {
            // KODU: Convert to checked return statement.
            return list.get(index);
        } catch (IndexOutOfBoundsException e) {
            return null;
        }
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
}

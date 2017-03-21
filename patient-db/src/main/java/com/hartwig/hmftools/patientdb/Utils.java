package com.hartwig.hmftools.patientdb;

import java.text.DateFormat;
import java.text.ParseException;
import java.util.Collection;
import java.util.Date;
import java.util.Iterator;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class Utils {
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
    static Date getDate(@Nullable String dateFieldValue, @NotNull DateFormat dateFormat) {
        try {
            return dateFormat.parse(dateFieldValue);
        } catch (ParseException | NullPointerException e) {
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
}

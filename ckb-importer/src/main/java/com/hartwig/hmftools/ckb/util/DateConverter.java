package com.hartwig.hmftools.ckb.util;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.Locale;

import org.jetbrains.annotations.Nullable;

public final class DateConverter {

    private static final DateTimeFormatter FORMAT = DateTimeFormatter.ofPattern("dd/MM/yyyy", Locale.ENGLISH);

    private DateConverter() {
    }

    @Nullable
    public static LocalDate toDate(@Nullable String string) {
        if (string == null) {
            return null;
        }

        try {
            return LocalDate.parse(string, FORMAT);
        } catch (DateTimeParseException e) {
            throw new IllegalStateException("Cannot convert string to date: " + string);
        }
    }
}

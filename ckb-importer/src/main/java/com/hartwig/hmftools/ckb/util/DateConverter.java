package com.hartwig.hmftools.ckb.util;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;

import org.jetbrains.annotations.Nullable;

public final class DateConverter {

    static final DateTimeFormatter FORMAT = DateTimeFormatter.ofPattern("MM/dd/yyyy");

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

package com.hartwig.hmftools.ckb.util;

import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Locale;

import org.jetbrains.annotations.Nullable;

public final class DateConverter {

    private static final DateFormat FORMAT = new SimpleDateFormat("dd/MM/yyyy", Locale.ENGLISH);

    private DateConverter() {
    }

    @Nullable
    public static Date toDate(@Nullable String string) {
        if (string == null) {
            return null;
        }

        try {
            return FORMAT.parse(string);
        } catch (ParseException e) {
            throw new IllegalStateException("Cannot convert string to date: " + string);
        }
    }
}

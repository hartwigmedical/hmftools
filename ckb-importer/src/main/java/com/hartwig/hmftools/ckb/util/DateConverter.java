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
    public static Date toDate(@Nullable String string) throws ParseException {
        return string != null ? FORMAT.parse(string) : null;
    }
}

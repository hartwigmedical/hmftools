package com.hartwig.hmftools.ckb.util;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Locale;

import org.jetbrains.annotations.Nullable;

import jdk.nashorn.internal.ir.annotations.Ignore;

public class DateConverter {

    private DateConverter() {

    }

    @Nullable
    @Ignore
    public static Date convertDate(@Nullable String dateString) throws ParseException {
        if (dateString == null) {
            return null;
        } else {
            return new SimpleDateFormat("dd/MM/yyyy", Locale.ENGLISH).parse(dateString);
        }
    }
}

package com.hartwig.hmftools.sullivan;

import org.jetbrains.annotations.NotNull;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

final class Logger {

    private static final DateFormat DTF = new SimpleDateFormat("hh:mm:ss");

    public static void log(@NotNull String info) {
        System.out.println(DTF.format(new Date()) + ": " + info);
    }
}

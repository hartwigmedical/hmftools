package com.hartwig.hmftools.ckb.dao;

import java.time.LocalDate;

import org.jetbrains.annotations.Nullable;

final class Util {

    private Util() {
    }

    @Nullable
    public static java.sql.Date sqlDate(@Nullable LocalDate date) {
        return date != null ? java.sql.Date.valueOf(date) : null;
    }

    public static byte toByte(boolean value) {
        return value ? (byte) 1 : (byte) 0;
    }
}

package com.hartwig.hmftools.iclusion.dao;

import org.jetbrains.annotations.NotNull;

public class Utils {

    public static final int DB_BATCH_INSERT_SIZE = 1000;

    private Utils() {
    }

    static boolean anyNull(@NotNull Object... arguments) {
        for (Object object : arguments) {
            if (object == null) {
                return true;
            }
        }
        return false;
    }
}

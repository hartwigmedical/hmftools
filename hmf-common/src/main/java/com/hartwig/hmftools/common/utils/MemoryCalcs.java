package com.hartwig.hmftools.common.utils;

import static java.lang.Math.round;

public final class MemoryCalcs
{
    private static final long MEGABYTE = 1024L * 1024L;

    public static int calcMemoryUsage()
    {
        Runtime runtime = Runtime.getRuntime();
        long memory = runtime.totalMemory() - runtime.freeMemory();
        return round(memory / MEGABYTE);
    }
}

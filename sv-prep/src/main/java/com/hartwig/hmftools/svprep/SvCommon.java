package com.hartwig.hmftools.svprep;

import static java.lang.Math.round;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvCommon
{
    public static final Logger SV_LOGGER = LogManager.getLogger(SvCommon.class);

    public static final String DELIM = ",";
    public static final String ITEM_DELIM = ";";
    public static final String SUB_ITEM_DELIM = ":";

    private static final long MEGABYTE = 1024L * 1024L;

    public static int calcMemoryUsage(boolean runGc)
    {
        Runtime runtime = Runtime.getRuntime();

        if(runGc)
            runtime.gc();

        long memory = runtime.totalMemory() - runtime.freeMemory();
        return round(memory / MEGABYTE);
    }
}

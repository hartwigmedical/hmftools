package com.hartwig.hmftools.sage;

import static java.lang.Math.round;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SageCommon
{
    public static final String DELIM = "\t";

    public static final Logger SG_LOGGER = LogManager.getLogger(SageCommon.class);

    private static final long MEGABYTE = 1024L * 1024L;

    public static int calcCurrentMemoryUsage(boolean runGc)
    {
        Runtime runtime = Runtime.getRuntime();

        if(runGc)
            runtime.gc();

        long memory = runtime.totalMemory() - runtime.freeMemory();
        return round(memory / MEGABYTE);
    }

}

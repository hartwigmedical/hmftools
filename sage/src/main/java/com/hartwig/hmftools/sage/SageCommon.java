package com.hartwig.hmftools.sage;

import com.hartwig.hmftools.common.utils.MemoryCalcs;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SageCommon
{
    public static final String DELIM = "\t";

    public static final Logger SG_LOGGER = LogManager.getLogger(SageCommon.class);

    public static int calcMemoryUsage(boolean runGc)
    {
        if(runGc)
            System.gc();

        return MemoryCalcs.calcMemoryUsage();
    }

    public static void logMemoryUsage(final SageConfig config, final String stage, int memory)
    {
        if(config.PerfWarnTime == 0)
            return;

        SG_LOGGER.info("stage({}) memory({})", stage, memory);
    }

}

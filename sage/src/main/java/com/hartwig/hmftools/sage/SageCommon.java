package com.hartwig.hmftools.sage;

import com.hartwig.hmftools.common.utils.MemoryCalcs;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SageCommon
{
    public static final String APP_NAME = "Sage";

    public static final String SAMPLE_DELIM = ",";

    public static final Logger SG_LOGGER = LogManager.getLogger(SageCommon.class);

    public static void logMemoryUsage(final double perfWarnTime, final String stage, int memory)
    {
        if(perfWarnTime == 0)
            return;

        SG_LOGGER.info("stage({}) memory({})", stage, memory);
    }

}

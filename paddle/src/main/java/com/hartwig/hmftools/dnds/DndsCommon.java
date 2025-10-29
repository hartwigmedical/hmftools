package com.hartwig.hmftools.dnds;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class DndsCommon
{
    public static final String APP_NAME = "DndsBuilder";

    public static final Logger DN_LOGGER = LogManager.getLogger(DndsCommon.class);

    public static final String SOMATIC_CACHE_DIR = "somatics";
    // public static final String SOMATIC_CACHE_DIR = "somatics";

    // constants
    public static final int MAX_REPEAT_COUNT = 7;
}

package com.hartwig.hmftools.dnds;

import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class DndsCommon
{
    public static final Logger DN_LOGGER = LogManager.getLogger(DndsCommon.class);

    public static final String SOMATIC_CACHE_DIR = "somatics";
    // public static final String SOMATIC_CACHE_DIR = "somatics";

    // constants
    public static final int MAX_REPEAT_COUNT = 7;

    public static void logVersion()
    {
        final VersionInfo version = new VersionInfo("dnds-builder.version");
        DN_LOGGER.info("DNDS Builder version: {}", version.version());

    }
}

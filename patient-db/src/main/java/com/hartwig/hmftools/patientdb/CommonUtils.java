package com.hartwig.hmftools.patientdb;

import com.hartwig.hmftools.common.utils.config.VersionInfo;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class CommonUtils
{
    public static final String APP_NAME = "PatientDb";

    public static final Logger LOGGER = LogManager.getLogger(CommonUtils.class);

    private CommonUtils() {}

    public static void logVersion()
    {
        final VersionInfo version = new VersionInfo("patientdb.version");
        LOGGER.info("Patient-DB version: {}", version.version());
    }

    static boolean anyNull(@NotNull Object... arguments)
    {
        for(Object object : arguments)
        {
            if(object == null)
            {
                return true;
            }
        }
        return false;
    }
}

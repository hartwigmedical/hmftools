package com.hartwig.hmftools.peach;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PeachUtils
{
    public static final Logger PCH_LOGGER = LogManager.getLogger(PeachApplication.class);

    public static final int GERMLINE_TOTAL_COPY_NUMBER = 2;
    public static final String COUNT_UNKNOWN_STRING = "UNKNOWN";

    @NotNull
    public static String convertCountToString(@Nullable Integer count)
    {
        return count == null ? COUNT_UNKNOWN_STRING : Integer.toString(count);
    }
}

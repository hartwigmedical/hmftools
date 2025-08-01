package com.hartwig.hmftools.wisp.common;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class CommonUtils
{
    public static final Logger CT_LOGGER = LogManager.getLogger(CommonUtils.class);

    public static final String APP_NAME = "Wisp";

    // common fields
    public static final String FLD_TUMOR_ID = "TumorId";
    public static final String FLD_CATEGORY = "Category";
    public static final String FLD_VARIANT = "Variant";

    public static final String BATCH_CONTROL_TAG = "batch_control";

    public static final int DEFAULT_PROBE_LENGTH = 120;
}

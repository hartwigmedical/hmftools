package com.hartwig.hmftools.purple;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class PurpleUtils
{
    public static final Logger PPL_LOGGER = LogManager.getLogger(PurpleUtils.class);

    public static String formatPurity(double purity) { return String.format("%.3f", purity); }

}

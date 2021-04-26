package com.hartwig.hmftools.purple;

import java.text.DecimalFormat;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class PurpleCommon
{
    public static final Logger PPL_LOGGER = LogManager.getLogger(PurpleCommon.class);

    public static final DecimalFormat formatDbl = new DecimalFormat("0.00");

    public static String formatPurity(double purity) { return String.format("%.3f", purity); }

}

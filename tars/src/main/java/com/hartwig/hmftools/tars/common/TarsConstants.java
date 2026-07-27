package com.hartwig.hmftools.tars.common;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class TarsConstants
{
    public static final Logger TARS_LOGGER = LogManager.getLogger(TarsConstants.class);

    public static final String APP_NAME = "Tars";

    // suffix appended to a chromosome to name its alt contig, e.g. chr1 -> chr1_tx
    public static final String ALT_CONTIG_SUFFIX = "_tx";
}

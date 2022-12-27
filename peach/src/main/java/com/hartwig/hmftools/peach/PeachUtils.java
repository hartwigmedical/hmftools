package com.hartwig.hmftools.peach;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class PeachUtils
{
    public static final Logger PCH_LOGGER = LogManager.getLogger(PeachApplication.class);

    public static final int GERMLINE_TOTAL_COPY_NUMBER = 2;

    public static final String TSV_DELIMITER = "\t";
    public static final String BED_FILE_DELIMITER = "\t";
    public static final String HAPLOTYPE_EVENT_DELIMITER = ",";
}

package com.hartwig.hmftools.common.rna;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class RnaCommon
{
    public static final String ISF_FILE_ID = ".isf.";

    // field names shared across various Isofox files
    public static final String FLD_GENE_SET_ID = "GeneSetId";
    public static final String FLD_FRAG_COUNT = "FragCount";
    public static final String FLD_DEPTH_START = "DepthStart";
    public static final String FLD_DEPTH_END = "DepthEnd";

    protected static final Logger RNA_LOGGER = LogManager.getLogger(RnaCommon.class);
}

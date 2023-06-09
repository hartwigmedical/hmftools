package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.utils.FileDelimiters.CSV_DELIM;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class RnaCommon
{
    public static final String ISF_FILE_ID = ".isf.";

    public static final String DELIMITER = CSV_DELIM;

    // field names shared across various Isofox files
    public static final String FLD_GENE_ID = "GeneId";
    public static final String FLD_GENE_SET_ID = "GeneSetId";
    public static final String FLD_GENE_NAME = "GeneName";
    public static final String FLD_CHROMOSOME = "Chromosome";
    public static final String FLD_TRANS_ID = "TransId";
    public static final String FLD_TRANS_NAME = "TransName";
    public static final String FLD_FRAG_COUNT = "FragCount";
    public static final String FLD_DEPTH_START = "DepthStart";
    public static final String FLD_DEPTH_END = "DepthEnd";

    protected static final Logger RNA_LOGGER = LogManager.getLogger(RnaCommon.class);
}

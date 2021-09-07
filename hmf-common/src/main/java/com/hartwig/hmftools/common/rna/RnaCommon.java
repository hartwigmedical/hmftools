package com.hartwig.hmftools.common.rna;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class RnaCommon
{
    public static final String ISF_FILE_ID = ".isf.";

    public static final String DELIMITER = ",";

    // field names shared across various Isofox files
    public static final String FLD_GENE_ID = "GeneId";
    public static final String FLD_GENE_SET_ID = "GeneSetId";
    public static final String FLD_GENE_NAME = "GeneName";
    public static final String FLD_CHROMOSOME = "Chromosome";
    public static final String FLD_TRANS_ID = "TransId";
    public static final String FLD_TRANS_NAME = "TransName";

    protected static final Logger RNA_LOGGER = LogManager.getLogger(RnaCommon.class);
}

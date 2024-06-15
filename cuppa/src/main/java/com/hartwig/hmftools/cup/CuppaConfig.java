package com.hartwig.hmftools.cup;

import com.hartwig.hmftools.common.utils.file.CommonFields;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CuppaConfig
{
    public static final String REF_DATA_DIR = "ref_data_dir";

    // reference data files
    public static final String REF_SAMPLE_DATA_FILE = "ref_sample_data_file";
    public static final String REF_SNV_COUNTS_FILE = "ref_snv_counts_file";
    public static final String REF_SNV_SAMPLE_POS_FREQ_FILE = "ref_sample_snv_pos_freq_file";
    public static final String REF_RNA_GENE_EXP_SAMPLE_FILE = "ref_gene_exp_sample_file";
    public static final String REF_RNA_ALT_SJ_SAMPLE_FILE = "ref_alt_sj_sample_file";

    public static final Logger CUP_LOGGER = LogManager.getLogger(CuppaConfig.class);

    // file fields
    public static final String FLD_SAMPLE_ID = CommonFields.FLD_SAMPLE_ID;
    public static final String FLD_CANCER_TYPE = CommonFields.FLD_CANCER_TYPE;
    public static final String FLD_CANCER_SUBTYPE = "CancerSubtype";
    public static final String FLD_RNA_READ_LENGTH = "RnaReadLength";
    public static final String CANCER_SUBTYPE_OTHER = "Other";
    public static final String DATA_DELIM = ",";
    public static final String SUBSET_DELIM = ";";
}

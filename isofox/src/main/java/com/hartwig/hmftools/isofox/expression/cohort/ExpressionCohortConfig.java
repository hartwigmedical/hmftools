package com.hartwig.hmftools.isofox.expression.cohort;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class ExpressionCohortConfig
{
    // cohort file generation and filters
    public final boolean WriteSampleGeneDistributionData;
    public final String CohortTransFile;
    public final String CancerTransFile;
    public final String GeneExpMatrixFile;

    public final double TpmThreshold;
    public final boolean UseLogTpm;
    public final double TpmRounding;

    public final String ExternalSource;
    public final boolean TranscriptScope;
    public final String CancerGeneFiles;
    public final boolean DistributionByCancerType;
    public final boolean LogElevatedDistributions;
    public final boolean ApplyTpmWriteLimit;

    private static final String EXTERNAL_SOURCE = "exp_external_source";
    private static final String EXTERNAL_COMPARE_TRANSCRIPTS = "exp_compare_transcripts";

    public static final String WRITE_SAMPLE_GENE_DISTRIBUTION_DATA = "write_sample_gene_dist";

    public static final String TPM_THRESHOLD = "tpm_threshold";
    public static final String USE_LOG_TPM = "use_log_tpm";
    public static final String TPM_ROUNDING = "tpm_rounding";
    public static final String APPLY_TPM_WRITE_THRESHOLD = "apply_tpm_write_threshold";

    public static final String GENE_EXP_MATRIX_FILE = "gene_exp_matrix_file";
    public static final String COHORT_TRANS_FILE = "cohort_trans_file";
    public static final String CANCER_TRANS_FILE = "cancer_trans_file";
    public static final String CANCER_GENE_FILES = "cancer_gene_files";
    public static final String DIST_BY_CANCER_TYPE = "dist_by_cancer";
    public static final String LOG_ELEVATED_DIST = "log_elevated_dist";

    public static final String SOURCE_ISOFOX = "ISOFOX";
    public static final String EXT_SOURCE_SALMON = "SALMON";
    public static final String EXT_SOURCE_RSEM = "RSEM";

    public ExpressionCohortConfig(final ConfigBuilder configBuilder)
    {
        WriteSampleGeneDistributionData = configBuilder.hasFlag(WRITE_SAMPLE_GENE_DISTRIBUTION_DATA);
        TpmThreshold = configBuilder.getDecimal(TPM_THRESHOLD);
        TpmRounding = configBuilder.getDecimal(TPM_ROUNDING);
        UseLogTpm = configBuilder.hasValue(USE_LOG_TPM);

        GeneExpMatrixFile = configBuilder.getValue(GENE_EXP_MATRIX_FILE);
        CohortTransFile = configBuilder.getValue(COHORT_TRANS_FILE);
        CancerTransFile = configBuilder.getValue(CANCER_TRANS_FILE);
        TranscriptScope = configBuilder.hasFlag(EXTERNAL_COMPARE_TRANSCRIPTS);

        ExternalSource = configBuilder.getValue(EXTERNAL_SOURCE);
        CancerGeneFiles = configBuilder.getValue(CANCER_GENE_FILES);
        DistributionByCancerType = configBuilder.hasFlag(DIST_BY_CANCER_TYPE);
        LogElevatedDistributions = configBuilder.hasFlag(LOG_ELEVATED_DIST);
        ApplyTpmWriteLimit = configBuilder.hasFlag(APPLY_TPM_WRITE_THRESHOLD);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(EXTERNAL_SOURCE, "List of sources to compare fusions between");
        configBuilder.addFlag(WRITE_SAMPLE_GENE_DISTRIBUTION_DATA, "Write per-sample gene distribution data file");
        configBuilder.addPath(GENE_EXP_MATRIX_FILE, false, "Cohort gene expression file");
        configBuilder.addPath(COHORT_TRANS_FILE, false, "Cohort transcript distribution file");
        configBuilder.addPath(CANCER_TRANS_FILE, false, "Cancer transcript distribution file");
        configBuilder.addConfigItem(EXTERNAL_COMPARE_TRANSCRIPTS, false, "Compare at transcript level, other default is by gene");
        configBuilder.addDecimal(TPM_ROUNDING, "TPM/FPM rounding factor, base-10 integer (2 = 1%)", 2);
        configBuilder.addDecimal(TPM_THRESHOLD, "Only write transcripts with TPM greater than this", 0);
        configBuilder.addFlag(USE_LOG_TPM, "Convert TPM to log");
        configBuilder.addConfigItem(CANCER_GENE_FILES, false, "Cancer gene distribution files, format: CancerType1-File1;CancerType2-File2");
        configBuilder.addFlag(DIST_BY_CANCER_TYPE, "Produce cancer gene distributions");
        configBuilder.addFlag(LOG_ELEVATED_DIST, "Log elevated gene distributions");
        configBuilder.addFlag(APPLY_TPM_WRITE_THRESHOLD, "Log elevated gene distributions");
    }
}

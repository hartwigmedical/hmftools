package com.hartwig.hmftools.isofox.expression.cohort;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

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

    public ExpressionCohortConfig(final CommandLine cmd)
    {
        WriteSampleGeneDistributionData = cmd.hasOption(WRITE_SAMPLE_GENE_DISTRIBUTION_DATA);
        TpmThreshold = Double.parseDouble(cmd.getOptionValue(TPM_THRESHOLD, "0"));
        TpmRounding = Double.parseDouble(cmd.getOptionValue(TPM_ROUNDING, "2"));
        UseLogTpm = cmd.hasOption(USE_LOG_TPM);

        GeneExpMatrixFile = cmd.getOptionValue(GENE_EXP_MATRIX_FILE);
        CohortTransFile = cmd.getOptionValue(COHORT_TRANS_FILE);
        CancerTransFile = cmd.getOptionValue(CANCER_TRANS_FILE);
        TranscriptScope = cmd.hasOption(EXTERNAL_COMPARE_TRANSCRIPTS);

        ExternalSource = cmd.getOptionValue(EXTERNAL_SOURCE);
        CancerGeneFiles = cmd.getOptionValue(CANCER_GENE_FILES);
        DistributionByCancerType = cmd.hasOption(DIST_BY_CANCER_TYPE);
        LogElevatedDistributions = cmd.hasOption(LOG_ELEVATED_DIST);
        ApplyTpmWriteLimit = cmd.hasOption(APPLY_TPM_WRITE_THRESHOLD);
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(EXTERNAL_SOURCE, true, "List of sources to compare fusions between");
        options.addOption(WRITE_SAMPLE_GENE_DISTRIBUTION_DATA, false, "Write per-sample gene distribution data file");
        options.addOption(GENE_EXP_MATRIX_FILE, true, "Cohort gene expression file");
        options.addOption(COHORT_TRANS_FILE, true, "Cohort transcript distribution file");
        options.addOption(CANCER_TRANS_FILE, true, "Cancer transcript distribution file");
        options.addOption(EXTERNAL_COMPARE_TRANSCRIPTS, false, "Compare at transcript level, other default is by gene");
        options.addOption(TPM_ROUNDING, true, "TPM/FPM rounding factor, base-10 integer (default=2, ie 1%)");
        options.addOption(TPM_THRESHOLD, true, "Only write transcripts with TPM greater than this");
        options.addOption(USE_LOG_TPM, false, "Convert TPM to log");
        options.addOption(CANCER_GENE_FILES, true, "Cancer gene distribution files, format: CancerType1-File1;CancerType2-File2");
        options.addOption(DIST_BY_CANCER_TYPE, false, "Produce cancer gene distributions");
        options.addOption(LOG_ELEVATED_DIST, false, "Log elevated gene distributions");
        options.addOption(APPLY_TPM_WRITE_THRESHOLD, false, "Log elevated gene distributions");
    }
}

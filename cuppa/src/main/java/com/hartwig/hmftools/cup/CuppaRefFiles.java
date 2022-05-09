package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SOMATIC_VCF_SUFFIX;
import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SV_VCF_SUFFIX;
import static com.hartwig.hmftools.cup.CuppaConfig.formSamplePath;

import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;

public class CuppaRefFiles
{
    public static final String CUP_REF_FILE_PREFIX = "cup_ref";

    private static String formatRefFilename(final String fileType)
    {
        return String.format("%s_%s.csv", CUP_REF_FILE_PREFIX, fileType);
    }

    public static final String REF_FILE_SAMPLE_DATA = formatRefFilename("sample_data");
    public static final String REF_FILE_SNV_COUNTS = formatRefFilename("snv_counts");
    public static final String REF_FILE_CANCER_POS_FREQ_COUNTS = formatRefFilename("cancer_pos_freq_counts");
    public static final String REF_FILE_SAMPLE_POS_FREQ_COUNTS = formatRefFilename("sample_pos_freq_counts");
    public static final String REF_FILE_CN_ADJUSTED_SAMPLE_POS_FREQ_COUNTS = formatRefFilename("sample_cn_adj_pos_freq");
    public static final String REF_FILE_COPY_NUMBER_PROFILE = formatRefFilename("copy_number_profile");
    public static final String REF_FILE_FEATURE_PREV = formatRefFilename("feature_prev");
    public static final String REF_FILE_DRIVER_AVG = formatRefFilename("driver_avg");
    public static final String REF_FILE_TRAIT_PERC = formatRefFilename("sample_trait_percentiles");
    public static final String REF_FILE_TRAIT_RATES = formatRefFilename("sample_trait_rates");
    public static final String REF_FILE_GENDER_RATES = formatRefFilename("gender_rates");
    public static final String REF_FILE_SIG_PERC = formatRefFilename("sig_percentiles");
    public static final String REF_FILE_SV_PERC = formatRefFilename("sv_percentiles");
    public static final String REF_FILE_GENE_EXP_CANCER = formatRefFilename("gene_exp_cancer");
    public static final String REF_FILE_GENE_EXP_SAMPLE = formatRefFilename("gene_exp_sample");
    public static final String REF_FILE_ALT_SJ_CANCER = formatRefFilename("alt_sj_cancer");
    public static final String REF_FILE_ALT_SJ_SAMPLE = formatRefFilename("alt_sj_sample");
    public static final String REF_FILE_SNV_SIGNATURES = formatRefFilename("snv_signatures");

    // cohort files for building reference data
    public static final String COHORT_REF_FILE_SV_DATA_FILE = formatRefFilename("cohort_sv_data");
    public static final String COHORT_REF_FILE_TRAITS_DATA_FILE = formatRefFilename("cohort_traits_data");
    public static final String COHORT_REF_FILE_SIG_DATA_FILE = formatRefFilename("cohort_signature_data");
    public static final String COHORT_REF_FILE_FEATURE_DATA_FILE = formatRefFilename("cohort_feature_data");

    public static final String purpleSomaticVcfFile(final String sampleDir, final String sampleId)
    {
        String basePath = formSamplePath(sampleDir, sampleId);
        return basePath + sampleId + PURPLE_SOMATIC_VCF_SUFFIX;
    }

    public static final String purpleSvFile(final String sampleDir, final String sampleId)
    {
        String basePath = formSamplePath(sampleDir, sampleId);
        return basePath + sampleId + PURPLE_SV_VCF_SUFFIX;
    }

    public static final String purplePurityFile(final String sampleDir, final String sampleId)
    {
        String basePath = formSamplePath(sampleDir, sampleId);
        return PurityContextFile.generateFilenameForReading(basePath, sampleId);
    }

    public static final String purpleCopyNumberFile(final String sampleDir, final String sampleId)
    {
        String basePath = formSamplePath(sampleDir, sampleId);
        return PurpleCopyNumberFile.generateFilenameForReading(basePath, sampleId);
    }
}

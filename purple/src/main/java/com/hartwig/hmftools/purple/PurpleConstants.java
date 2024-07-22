package com.hartwig.hmftools.purple;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class PurpleConstants
{
    // constants prefixed with 'DEFAULT' can be overridden in config

    // common
    public static final int WINDOW_SIZE = 1000;
    public static final double BAF_PNT_5 = 0.5;

    // TMB calcs
    public static final double MB_PER_GENOME = 2859;
    public static final double CODING_BASES_PER_GENOME = 3.188e7; // calculated from GRCh38 canonical transcripts (overlaps ignored)

    // no-tumor
    public static final int NO_TUMOR_BAF_TOTAL = 3000;
    public static final double NO_TUMOR_DEPTH_RATIO_MIN = 0.8;
    public static final double NO_TUMOR_DEPTH_RATIO_MAX = 1.2;

    // target regions
    public static final int TARGET_REGIONS_MAX_DELETED_GENES = 1500;
    public static final double ASSUMED_BIALLELIC_FRACTION = 0.2;

    // purity fitting
    public static final double MAX_DIPLOID_COPY_NUMBER = 1.2;
    public static final double MIN_DIPLOID_COPY_NUMBER = 0.8;

    // segmentation and regions
    public static final int CENTROMERIC_WIDTH = 4_000_000;
    public static final double GERMLINE_AMP_RATIO = 1.3;
    public static final double GERMLINE_DEL_RATIO = 0.7;
    public static final int GERMLINE_DEL_MIN_LENGTH = 5_000_000;
    public static final int GERMLINE_AMP_DEL_EXCLUSION_CHR_1 = 50_000_000;
    public static final int GERMLINE_AMP_DEL_EXCLUSION_CHR_9 = 135_000_000;
    public static final int GERMLINE_AMP_DEL_EXCLUSION_CHR_17 = 75_000_000;
    public static final int GERMLINE_AMP_DEL_EXCLUSION_CHR_19 = 20_000_000;

    // purity calcs
    public static final double MIN_PURITY_DEFAULT = 0.08;
    public static final double MAX_PURITY_DEFAULT = 1.0;
    public static final double PURITY_INCREMENT_DEFAULT = 0.01;
    public static final double MIN_PLOIDY_DEFAULT = 1.0;
    public static final double MAX_PLOIDY_DEFAULT = 8;

    public static final int MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT = 30;
    public static final int MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT = 150;
    public static final int TARGETED_MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT = 3;

    public static final double PLOIDY_PENALTY_FACTOR_DEFAULT = 0.4;

    public static final double PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT = 0.05;
    public static final double TARGETED_PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT = 0.1;

    public static final double PLOIDY_PENALTY_MIN_STANDARD_DEVIATION_DEFAULT = 1.5;

    public static final double PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT = 1;
    public static final double TARGETED_PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT = 3;

    public static final double PLOIDY_PENALTY_SUB_MIN_ADDITIONAL_DEFAULT = 1.5;

    public static final double PLOIDY_PENALTY_MIN_DEFAULT = 0.1;
    public static final double TARGETED_PLOIDY_PENALTY_MIN_DEFAULT = 0.2;

    public static final double TARGETED_DEVIATION_PENALTY_GC_MIN_ADJUST_DEFAULT = 0.25;
    public static final double TARGETED_GC_RATIO_EXPONENT_DEFAULT = 3;

    // somatic fitting
    public static final double SNV_HOTSPOT_VAF_PROBABILITY = 0.01;
    public static final int SNV_HOTSPOT_MAX_SNV_COUNT = 2000;
    public static final double SNV_FITTING_MAPPABILITY = 1.0;
    public static final int SNV_FITTING_MAX_REPEATS = 3;
    public static final double SOMATIC_MIN_PURITY_DEFAULT = 0.17;
    public static final double SOMATIC_MIN_PURITY_SPREAD_DEFAULT = 0.15;
    public static final int SOMATIC_MIN_PEAK_DEFAULT = 4;
    public static final int SOMATIC_MIN_VARIANTS_DEFAULT = 10;
    public static final double SOMATIC_PENALTY_WEIGHT_DEFAULT = 1.5;
    public static final double HIGHLY_DIPLOID_PERCENTAGE_DEFAULT = 0.97;

    // tumor-only somatic fitting
    public static final double SOMATIC_FIT_TUMOR_ONLY_PURITY_MIN = 0.92;
    public static final double SOMATIC_FIT_TUMOR_ONLY_PLOIDY_MIN = 1.8;
    public static final double SOMATIC_FIT_TUMOR_ONLY_PLOIDY_MAX = 2.2;
    public static final double SOMATIC_FIT_TUMOR_ONLY_MIN_VAF = 0.04;

    // somatic fitting readjustment
    public static final double SNV_READJUST_CN_MIN = 1.8;
    public static final double SNV_READJUST_CN_MAX = 2.2;
    public static final double SNV_READJUST_MINOR_ALLELE_MIN = 0.5;
    public static final double SNV_READJUST_PROB_THRESHOLD = 0.005;
    public static final int SNV_READJUST_PROB_THRESHOLD_MIN_VARIANTS = 5;
    public static final double SNV_READJUST_EXPECTED_VARIANT_RATIO = 4;
    public static final double SNV_READJUST_EXPECTED_VARIANT_COUNT_LOW_PURITY_LEVEL = 0.2;
    public static final int SNV_READJUST_EXPECTED_VARIANT_COUNT_LOW_PURITY = 10;
    public static final int SNV_READJUST_EXPECTED_VARIANT_COUNT = 20;
    public static final double SNV_READJUST_PURITY_INCREMENT = 0.005;
    public static final double INVALID_PURITY = -1;

    public static final double CLONALITY_BIN_WIDTH = 0.05;
    public static final double CLONALITY_MAX_PLOIDY = 10;
    public static final int MIN_TOTAL_SV_FRAGMENT_COUNT = 1000;
    public static final int MIN_TOTAL_SOMATIC_VAR_ALLELE_READ_COUNT = 1000;

    // SV recovery
    public static final int DEFAULT_RECOVERY_MIN_MATE_QUAL_SCORE = 300;
    public static final int DEFAULT_RECOVERY_MIN_SGL_QUAL_SCORE = 500;
    public static final double RECOVERY_MIN_LENGTH = 1000;
    public static final double RECOVERY_MIN_PLOIDY = 0.5;
    public static final double RECOVERY_MIN_PLOIDY_PERC = 0.5;
    public static final int RECOVERY_MIN_MATE_UNCERTAINTY = 150;
    public static final int RECOVERY_UNBALANCED_MIN_DEPTH_WINDOW_COUNT = 5;
    public static final double RECOVERY_UNBALANCED_MIN_UNEXPLAINED_COPY_NUMBER_CHANGE = 0.6;
    public static final double RECOVERY_UNBALANCED_MIN_UNEXPLAINED_COPY_NUMBER_CHANGE_PERC = 0.2;

    // germline deletions
    public static final int GERMLINE_DEL_GENE_BUFFER = 500;
    public static final int GERMLINE_DEL_REGION_MIN = 1000;
    public static final int GERMLINE_DEL_REGION_MATCH_BUFFER = 2000;
    public static final double GERMLINE_DEL_CN_CONSISTENCY_MIN = 0.5;
    public static final double GERMLINE_DEL_CN_CONSISTENCY_MACN_PERC = 0.2;
    public static final double GERMLINE_DEL_NORMAL_RATIO = 0.65;
    public static final int GERMLINE_DEL_COHORT_FREQ = 4;

    // copy number smoothing
    public static final double MIN_OBSERVED_BAF_CHANGE = 0.03;
    public static final double MAX_DEVIATION_ADJUSTMENT = 0.20;
    public static final double MIN_ABSOLUTE_COPY_NUMBER_TOLERANCE = 0.3;
    public static final double MIN_RELATIVE_COPY_NUMBER_TOLERANCE = 0.12;
    public static final double MIN_ABSOLUTE_COPY_NUMBER_ADDITION = 1;
    public static final double MIN_RELATIVE_COPY_NUMBER_ADDITION = 0.8;

    public static final ChrBaseRegion CDKN2A_DELETION_REGION = new ChrBaseRegion("9", 9000000, 12000000);
    public static final double MAX_SOMATIC_FIT_DELETED_PERC = 0.003;

    // somatic subclonality peaks
    public static final int PEAK_BIN_COUNT = 10;
    public static final double PEAK_BIN_WIDTH = 0.01;
    public static final double PEAK_BIN_MIN_AVERAGE_WEIGHT = 0.4;
    public static final double PEAK_BIN_CLONAL_PLOIDY = 0.85;
    public static final double MAX_UNEXPLAINED_WEIGHT_PERCENT = 0.01;

    // drivers
    public static final int MAX_INDEL_DRIVER_REPEAT_COUNT = 7;

    // biallelic (no wildtype allele remaining) probability
    public static final int BIALLELIC_LOH_GROWTH_RATE = 40;
    public static final double BIALLELIC_THRESHOLD_PARAMETER_I = 1.4;
    public static final double BIALLELIC_THRESHOLD_PARAMETER_II = 0.8;
    public static final double BIALLELIC_LOH_BASE_ERROR_RATE = 0.02;
}

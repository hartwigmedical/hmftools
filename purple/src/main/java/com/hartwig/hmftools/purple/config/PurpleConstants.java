package com.hartwig.hmftools.purple.config;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class PurpleConstants
{
    // constants prefixed with 'DEFAULT' can be overridden in config

    // common
    public static final int WINDOW_SIZE = 1000;

    // TMB calcs
    public static final double MB_PER_GENOME = 2859;
    public static final double CODING_BASES_PER_GENOME = 3.188e7; // calculated from GRCh38 canonical transcripts (overlaps ignored)

    // no-tumor
    public static final int NO_TUMOR_BAF_TOTAL = 3000;
    public static final double NO_TUMOR_DEPTH_RATIO_MIN = 0.8;
    public static final double NO_TUMOR_DEPTH_RATIO_MAX = 1.2;

    public static final int TARGET_REGIONS_MAX_DELETED_GENES = 1500;

    // somatic fitting
    public static final double SNV_HOTSPOT_VAF_PROBABILITY = 0.01;
    public static final int SNV_HOTSPOT_MAX_SNV_COUNT = 2000;
    public static final double SNV_FITTING_MAPPABILITY = 1.0;
    public static final int SNV_FITTING_MAX_REPEATS = 3;
    public static final double SOMATIC_MIN_PURITY_DEFAULT = 0.17;
    public static final double SOMATIC_MIN_PURITY_SPREAD_DEFAULT = 0.15;
    public static final int SOMATIC_MIN_PEAK_DEFAULT = 4;
    public static final int SOMATIC_MIN_VARIANTS_DEFAULT = 10;
    public static final double SOMATIC_PENALTY_WEIGHT_DEFAULT = 1;
    public static final double HIGHLY_DIPLOID_PERCENTAGE_DEFAULT = 0.97;

    public static final double CLONALITY_BIN_WIDTH = 0.05;
    public static final double CLONALITY_MAX_PLOIDY = 10;
    public static final int MIN_TOTAL_SV_FRAGMENT_COUNT = 1000;
    public static final int MIN_TOTAL_SOMATIC_VAR_ALLELE_READ_COUNT = 2000;


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


}

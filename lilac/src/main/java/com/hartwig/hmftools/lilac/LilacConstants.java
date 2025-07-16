package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.lilac.GeneCache.longGeneName;

import java.util.List;

import com.google.common.collect.Lists;

public final class LilacConstants
{
    private LilacConstants() {}

    public static final String APP_NAME = "Lilac";

    public static final byte DEFAULT_MIN_BASE_QUAL = 30;
    public static byte LOW_BASE_QUAL_THRESHOLD = DEFAULT_MIN_BASE_QUAL; // may be adjusted by config or dynamically from median quals

    public static final double LOW_BASE_TRIM_PERC = 0.35;
    public static final double MAX_LOW_BASE_PERC = 0.5;
    public static final int DEFAULT_FRAGS_PER_ALLELE = 7;
    public static final int DEFAULT_FRAGS_REMOVE_SGL = 40;
    public static final double DEFAULT_TOP_SCORE_THRESHOLD = 0.005;

    public static final double DEFAULT_MIN_EVIDENCE_FACTOR = 0.006;
    public static final double DEFAULT_MIN_HIGH_QUAL_EVIDENCE_FACTOR = 0.003;
    public static final int DEFAULT_MIN_EVIDENCE_SUPPORT = 1;
    public static final int DEFAULT_MIN_DEPTH_FILTER = 10;

    // can be overridden in config
    public static double MIN_EVIDENCE_FACTOR = DEFAULT_MIN_EVIDENCE_FACTOR;
    public static double MIN_HIGH_QUAL_EVIDENCE_FACTOR = DEFAULT_MIN_HIGH_QUAL_EVIDENCE_FACTOR;
    public static int MIN_EVIDENCE_SUPPORT = DEFAULT_MIN_EVIDENCE_SUPPORT;
    public static int MIN_DEPTH_FILTER = DEFAULT_MIN_DEPTH_FILTER;

    public static final int DEFAULT_MAX_REF_FRAGMENTS = 1000;

    public static final double BASE_QUAL_PERCENTILE = 0.25;

    public static final double DEFAULT_HLA_Y_FRAGMENT_THRESHOLD = 0.005;
    public static final double COMMON_ALLELES_FREQ_CUTOFF = 0.001;
    public static final double RARE_ALLELES_FREQ_CUTOFF = 0.0001;

    // candidate filtering
    public static final double MIN_HI_CONF_UNIQUE_GROUP_COVERAGE = 0.02;
    public static final double MIN_HI_CONF_UNIQUE_PROTEIN_COVERAGE = 0.01;
    public static final double MIN_LOW_CONF_GROUP_UNIQUE_COVERAGE = 0.001;
    public static final double MIN_LOW_CONF_GROUP_TOTAL_COVERAGE = 0.035;
    public static final int MIN_WILDCARD_FRAGMENTS = 2;

    // scoring of allele combinations
    public static final double MIN_POPULATION_FREQUENCY = 0.0001;
    public static final double FREQUENCY_SCORE_PENALTY = 0.0018;
    public static final double SOLUTION_COMPLEXITY_PENALTY_WEIGHT = 0.002;
    public static final double RECOVERY_SCORE_PENALTY = 0;

    // QC thresholds
    public static int FAIL_LOW_COVERAGE_THRESHOLD = 360;
    public static double WARN_LOW_COVERAGE_THRESHOLD = 180;
    public static final double WARN_LOW_BASE_QUAL_THRESHOLD = 25;
    public static final double WARN_UNMATCHED_HAPLOTYPE_SUPPORT = 0.01;
    public static final int LOG_UNMATCHED_HAPLOTYPE_SUPPORT = 3;
    public static final double WARN_INDEL_THRESHOLD = 0.005;
    public static final double WARN_LOW_COVERAGE_DEPTH = 10;
    public static final int DEFAULT_FATAL_TOTAL_LOW_COVERAGE_POSITIONS = 300;

    public static final String HLA_PREFIX = "HLA-";
    public static String HLA_CHR = "6"; // note this is set as a versioned chromosome during initialisation

    // MHC class 1 constants
    public static final String GENE_A = "A";
    public static final String GENE_B = "B";
    public static final String GENE_C = "C";
    public static final String GENE_Y = "Y";
    public static final String GENE_H = "H";

    public static final String HLA_A = longGeneName(GENE_A);
    public static final String HLA_B = longGeneName(GENE_B);
    public static final String HLA_C = longGeneName(GENE_C);

    public static final List<String> CLASS_1_EXCLUDED_ALLELES = Lists.newArrayList(
            "A*31:135", "A*33:191", "A*02:783", "B*07:282",

            // Similar to HLA-Y
            "A*30:205", "A*30:207", "A*30:225", "A*30:228", "A*01:81", "A*01:237");

    // common INDEL associated with allele C*04:09N
    public static final String STOP_LOSS_ON_C_ALLELE = "C*04:09N";

    public static final int SPLICE_VARIANT_BUFFER = 5;

    public static final int COMPLEX_PERMS_THRESHOLD = 100000;

    // output file IDs
    public static final String LILAC_FILE_ID = ".lilac.";

    public static final String LILAC_FILE_CANDIDATE_COVERAGE = "candidates.coverage.tsv";
    public static final String LILAC_FILE_CANDIDATE_FRAGS = "candidates.fragments.tsv";
    public static final String LILAC_FILE_FRAGMENTS = "fragments.tsv";
    public static final String LILAC_FILE_READS = "reads.tsv";

    public static final String LILAC_FILE_CANDIDATE_AA = "candidates.aminoacids.txt";
    public static final String LILAC_FILE_CANDIDATE_NUC = "candidates.nucleotides.txt";

    public static final String LILAC_FILE_SOMATIC_VCF = "somatic.vcf.gz";
    public static final String LILAC_FILE_HLA_Y_COVERAGE = "hlay.coverage.tsv";
    public static final String LILAC_FILE_HLA_Y_FRAGMENTS = "hlay.fragments.tsv";
}

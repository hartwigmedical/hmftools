package com.hartwig.hmftools.lilac;

import java.util.List;
import java.util.OptionalDouble;

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

    public static final double DEFAULT_HLA_Y_FRAGMENT_THRESHOLD = 0.003;
    public static final double COMMON_ALLELES_FREQ_CUTOFF = 0.005;
    public static final double RARE_ALLELES_FREQ_CUTOFF = 0.0001;

    // candidate filtering
    public static final double MIN_HI_CONF_UNIQUE_GROUP_COVERAGE = 0.02;
    public static final double MIN_HI_CONF_UNIQUE_PROTEIN_COVERAGE = 0.01;
    public static final double MIN_LOW_CONF_GROUP_UNIQUE_COVERAGE = 0.001;
    public static final double MIN_LOW_CONF_GROUP_TOTAL_COVERAGE = 0.035;
    public static final int MIN_WILDCARD_FRAGMENTS = 2;

    // scoring of allele combinations
    public static final double MIN_POPULATION_FREQUENCY = 0.0001;

    public static final double DEFAULT_FREQUENCY_SCORE_PENALTY = 0.0009;
    public static double FREQUENCY_SCORE_PENALTY = DEFAULT_FREQUENCY_SCORE_PENALTY;

    public static OptionalDouble SOLUTION_COMPLEXITY_PENALTY_WEIGHT = OptionalDouble.empty();
    public static final int EXON_CHUNK_SIZE = 46;
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

    public static String HLA_CHR = "6"; // note this is set as a versioned chromosome during initialisation
    public static final String HLA_PREFIX = "HLA-";

    public static final List<String> CLASS_1_EXCLUDED_ALLELES = Lists.newArrayList(
            "A*31:135", "A*33:191", "A*02:783", "B*07:282",

            // Similar to HLA-Y
            "A*30:205", "A*30:207", "A*30:225", "A*30:228", "A*01:81", "A*01:237",

            // Similar to HLA-H
            "B*40:278", "A*03:423");

    public static final List<String> HLA_DRB1_EXCLUDED_ALLELES = Lists.newArrayList(
            "DRB1*14:242", // Similar to Exon 3 of DRB3
            "DRB1*14:141", // Similar to Exon 2 of DRB3
            "DRB1*14:262", // Similar to Exon 2 (1st half) of DRB3
            "DRB1*12:57" // Similar to Exon 2 (2nd half) of DRB3
    );

    // common INDEL associated with allele C*04:09N
    // TODO(mkcmkc): If we update the resources change to "C*04:09L"
    public static final String STOP_LOSS_ON_C_ALLELE = "C*04:09N";

    public static final int SPLICE_VARIANT_BUFFER = 5;

    public static final int COMPLEX_PERMS_THRESHOLD = 100000;

    // output file IDs
    public static final String LILAC_FILE_ID = ".lilac.";

    public static final String LILAC_FILE_CANDIDATE_COVERAGE = "candidates.coverage.tsv";
    public static final String LILAC_FILE_CANDIDATE_FRAGS = "candidates.fragments.tsv";
    public static final String LILAC_FILE_FRAGMENTS = "fragments.tsv";

    public static final String LILAC_FILE_CANDIDATE_AA = "candidates.aminoacids.txt";
    public static final String LILAC_FILE_CANDIDATE_NUC = "candidates.nucleotides.txt";

    public static final String LILAC_FILE_SOMATIC_VCF = "somatic.vcf.gz";
    public static final String LILAC_FILE_HLA_Y_COVERAGE = "HLA_Y.coverage.tsv";
    public static final String LILAC_FILE_HLA_Y_FRAGMENTS = "HLA_Y.fragments.tsv";
}

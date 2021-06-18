package com.hartwig.hmftools.lilac;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.read.Indel;

public class LilacConstants
{
    public static final int DEFAULT_MIN_BASE_QUAL = 30;
    public static final int DEFAULT_MIN_EVIDENCE = 2;
    public static final int DEFAULT_FRAGS_PER_ALLELE = 7;
    public static final int DEFAULT_FRAGS_REMOVE_SGL = 40;
    public static final double DEFAULT_TOP_SCORE_THRESHOLD = 0.005;
    public static final double DEFAULT_MIN_HIGH_QUAL_EVIDENCE_FACTOR = 0.000375;
    public static final double DEFAULT_MIN_EVIDENCE_FACTOR = 0.00075;

    // values applied as percentages
    public static final double COMMON_ALLELES_FREQ_CUTOFF = 0.001;
    public static final double MIN_CONF_UNIQUE_GROUP_COVERAGE = 0.02;
    public static final double MIN_CONF_UNIQUE_PROTEIN_COVERAGE = 0.01;
    public static final double HLA_Y_FRAGMENT_THRESHOLD = 0.01;
    public static final int MIN_WILDCARD_FRAGMENTS = 2;

    public static final double FREQUENCY_SCORE_PENALTY = 0.0018;
    public static final double HOMOZYGOUS_SCORE_PENALTY = 0.0036;
    public static final double RECOVERY_SCORE_PENALTY = 0.0055;
    public static final double WILDCARD_SCORE_PENALTY = 0.000015;

    // warning thresholds
    public static int FAIL_LOW_COVERAGE_THRESHOLD = 200;
    public static double WARN_LOW_COVERAGE_THRESHOLD = 50;
    public static final double WARN_LOW_BASE_QUAL_THRESHOLD = 25;
    public static final double WARN_UNMATCHED_HAPLOTYPE_SUPPORT = 0.01;
    public static final int LOG_UNMATCHED_HAPLOTYPE_SUPPORT = 3;
    public static final double WARN_INDEL_THRESHOLD = 0.005;
    public static final double WARN_LOW_COVERAGE_DEPTH = 10;
    public static final double FATAL_LOW_COVERAGE_THRESHOLD = 300;

    public static final String GENE_A = "A";
    public static final String GENE_B = "B";
    public static final String GENE_C = "C";
    public static final String GENE_Y = "Y";
    public static final String GENE_H = "H";

    public static final String HLA_PREFIX = "HLA-";
    public static final String HLA_A = longGeneName(GENE_A);
    public static final String HLA_B = longGeneName(GENE_B);
    public static final String HLA_C = longGeneName(GENE_C);

    public static final List<String> GENE_IDS = Lists.newArrayList(GENE_A, GENE_B, GENE_C);
    public static final List<String> HLA_GENES = Lists.newArrayList(HLA_A, HLA_B, HLA_C);
    public static String HLA_CHR = "6";

    public static final int EXPECTED_ALLELE_COUNT = 6;

    public static final List<String> EXCLUDED_ALLELES = Lists.newArrayList("A*31:135", "A*33:191");
    // previous PON list
    // A*01:81", "A*01:237", "A*11:126", "A*11:353", "A*25:68", "A*30:95", "A*30:136", "A*31:135", "A*33:191"

    // common INDEL associated with allele C*04:09N
    public static final String STOP_LOSS_ON_C_ALLELE = "C*04:09N";

    public static final List<Integer> A_EXON_BOUNDARIES = Lists.newArrayList(24, 114, 206, 298, 337, 348, 364);
    public static final List<Integer> B_EXON_BOUNDARIES = Lists.newArrayList(24, 114, 206, 298, 337, 348);
    public static final List<Integer> C_EXON_BOUNDARIES = Lists.newArrayList(24, 114, 206, 298, 338, 349, 365);

    public static final int NUC_LENGTH_A = 1098;
    public static final int NUC_LENGTH_B = 1089;
    public static final int NUC_LENGTH_C = 1101;

    public static final int MAX_AMINO_ACID_BOUNDARY = 298;

    public static final Map<String,List<Integer>> NUCLEOTIDE_EXON_BOUNDARIES = Maps.newHashMap();

    public static final int COMPLEX_PERMS_THRESHOLD = 100000;

    // technical
    public static final String DELIM = ",";
    public static final String ITEM_DELIM = ";";

    // common routines using constants
    public static List<Integer> getAminoAcidExonBoundaries(final String gene)
    {
        return gene.equals(GENE_A) ? A_EXON_BOUNDARIES : (gene.equals(GENE_B) ? B_EXON_BOUNDARIES : C_EXON_BOUNDARIES);
    }

    public static List<Integer> getNucleotideExonBoundaries(final String gene)
    {
        return NUCLEOTIDE_EXON_BOUNDARIES.get(gene);
    }

    public static String shortGeneName(final String gene)
    {
        return gene.startsWith(HLA_PREFIX) ? gene.substring(gene.length() - 1) : gene;
    }

    public static String longGeneName(final String gene)
    {
        return gene.length() == 1 ? HLA_PREFIX + gene : gene;
    }
}

package com.hartwig.hmftools.lilac;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.read.Indel;

public class LilacConstants
{
    public static final int DEFAULT_MIN_BASE_QUAL = 30;
    public static final int DEFAULT_MIN_EVIDENCE = 2;
    public static final int DEFAULT_FRAGS_PER_ALLELE = 7;
    public static final int DEFAULT_FRAGS_REMOVE_SGL = 40;
    public static final double DEFAULT_TOP_SCORE_THRESHOLD = 0.005;

    // values applied as percentages
    public static final double COMMON_ALLELES_FREQ_CUTOFF = 0.01;
    public static final double MIN_CONF_UNIQUE_GROUP_COVERAGE = 0.02;
    public static final double MIN_CONF_UNIQUE_PROTEIN_COVERAGE = 0.01;
    public static final double HLA_Y_FRAGMENT_THRESHOLD = 0.01;

    public static final double FREQUENCY_SCORE_PENALTY = 0.0015;
    public static final double HOMOZYGOUS_SCORE_PENALTY = 0.0045;
    public static final double RECOVERY_SCORE_PENALTY = 0.005;

    // warning thresholds
    public static final double WARN_UNMATCHED_HAPLOTYPE_SUPPORT = 0.008;
    public static final int LOG_UNMATCHED_HAPLOTYPE_SUPPORT = 3;
    public static final double WARN_INDEL_THRESHOLD = 0.004;

    public static final String GENE_A = "A";
    public static final String GENE_B = "B";
    public static final String GENE_C = "C";
    public static final String GENE_Y = "Y";

    public static final String HLA_PREFIX = "HLA-";
    public static final String HLA_A = longGeneName(GENE_A);
    public static final String HLA_B = longGeneName(GENE_B);
    public static final String HLA_C = longGeneName(GENE_C);

    public static final List<String> GENE_IDS = Lists.newArrayList(GENE_A, GENE_B, GENE_C);
    public static final List<String> HLA_GENES = Lists.newArrayList(HLA_A, HLA_B, HLA_C);
    public static final String HLA_CHR = "6";

    public static final int EXPECTED_ALLELE_COUNT = 6;

    public static final List<String> EXCLUDED_ALLELES = Lists.newArrayList(
            "A*01:81", "A*01:237", "A*11:353", "A*25:68", "A*30:95", "A*30:136", "A*31:135", "A*33:191");

    // common INDEL associated with allele C*04:09N
    public static final Indel STOP_LOSS_ON_C = new Indel("6", 31237115, "CN", "C");

    public static final List<Integer> A_EXON_BOUNDARIES = Lists.newArrayList(24, 114, 206, 298, 337, 348, 364);
    public static final List<Integer> B_EXON_BOUNDARIES = Lists.newArrayList(24, 114, 206, 298, 337, 348);
    public static final List<Integer> C_EXON_BOUNDARIES = Lists.newArrayList(24, 114, 206, 298, 338, 349, 365);

    public static final int MAX_AMINO_ACID_BOUNDARY = 298;

    public static final List<Integer> ALL_PROTEIN_EXON_BOUNDARIES = Lists.newArrayList(
            24, 114, 206, 298, 337, 348, 364, 338, 349, 365);

    public static final int COMPLEX_PERMS_THRESHOLD = 100000;

    // technical
    public static final String DELIM = ",";
    public static final String ITEM_DELIM = ";";

    public static List<Integer> getExonGeneBoundries(final String gene)
    {
        if(gene.equals(GENE_A))
            return A_EXON_BOUNDARIES;
        else if(gene.equals(GENE_B))
            return B_EXON_BOUNDARIES;
        else
            return C_EXON_BOUNDARIES;
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

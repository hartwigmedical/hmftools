package com.hartwig.hmftools.lilac;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class LilacConstants
{
    public static final int DEFAULT_MIN_BASE_QUAL = 30;
    public static final int DEFAULT_MIN_EVIDENCE = 2;
    public static final int DEFAULT_FRAGS_PER_ALLELE = 6;
    public static final int DEFAULT_FRAGS_REMOVE_SGL = 40;
    public static final int DEFAULT_MIN_CONF_UNIQUE_COVERAGE = 10;
    public static final int DEFAULT_MAX_DIST_FROM_TOP_SCORE = 3;

    public static final String HLA_A = "HLA-A";
    public static final String HLA_B = "HLA-B";
    public static final String HLA_C = "HLA-C";

    public static final Set<String> HLA_GENES = Sets.newHashSet(HLA_A, HLA_B, HLA_C);

    public static final HlaAllele DEFLATE_TEMPLATE = HlaAllele.fromString("A*01:01");

    public static final List<String> EXCLUDED_ALLELES = Lists.newArrayList(
            "A*01:81", "A*01:237", "A*33:191", "A*11:353", "A*30:95", "A*30:136", "A*31:135");

    public static final List<Integer> A_EXON_BOUNDARIES = Lists.newArrayList(24, 114, 206, 298, 337, 348, 364);
    public static final List<Integer> B_EXON_BOUNDARIES = Lists.newArrayList(24, 114, 206, 298, 337, 348);
    public static final List<Integer> C_EXON_BOUNDARIES = Lists.newArrayList(24, 114, 206, 298, 338, 349, 365);

    public static final int MAX_AMINO_ACID_BOUNDARY = 298;

    public static final List<Integer> ALL_PROTEIN_EXON_BOUNDARIES = Lists.newArrayList(
            24, 114, 206, 298, 337, 348, 364, 338, 349, 365);

    public static final List<Integer> ALL_NUCLEOTIDE_EXON_BOUNDARIES = Lists.newArrayList();

    public static final List<HmfTranscriptRegion> HLA_TRANSCRIPTS = Lists.newArrayList();

    public static final LociPosition LOCI_POSITION = new LociPosition();

    public static final String WILD_STR = "*";
    public static final char WILD_CHAR = WILD_STR.charAt(0);
    public static final String DEL_STR = ".";
    public static final char DEL_CHAR = DEL_STR.charAt(0);
}

package com.hartwig.hmftools.lilac;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.LociPosition;

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

    public static final com.hartwig.hmftools.lilac.hla.HlaAllele DEFLATE_TEMPLATE = com.hartwig.hmftools.lilac.hla.HlaAllele.from("A*01:01");

    public static final Set<com.hartwig.hmftools.lilac.hla.HlaAllele> EXCLUDED_ALLELES = Sets.newHashSet(
            com.hartwig.hmftools.lilac.hla.HlaAllele.from("A*01:81"),
            com.hartwig.hmftools.lilac.hla.HlaAllele.from("A*01:237"),
            com.hartwig.hmftools.lilac.hla.HlaAllele.from("A*33:191"),
            com.hartwig.hmftools.lilac.hla.HlaAllele.from("A*11:353"),
            com.hartwig.hmftools.lilac.hla.HlaAllele.from("A*30:95"),
            com.hartwig.hmftools.lilac.hla.HlaAllele.from("A*30:136"),
            HlaAllele.from("A*31:135"));

    public static final Set<Integer> A_EXON_BOUNDARIES = Sets.newHashSet(24, 114, 206, 298, 337, 348, 364);
    public static final Set<Integer> B_EXON_BOUNDARIES = Sets.newHashSet(24, 114, 206, 298, 337, 348);
    public static final Set<Integer> C_EXON_BOUNDARIES = Sets.newHashSet(24, 114, 206, 298, 338, 349, 365);

    public static final int MAX_AMINO_ACID_BOUNDARY = 298;

    public static final Set<Integer> ALL_PROTEIN_EXON_BOUNDARIES = Sets.newHashSet(
            24, 114, 206, 298, 337, 348, 364, 338, 349, 365);

    public static final List<Integer> ALL_NUCLEOTIDE_EXON_BOUNDARIES = Lists.newArrayList();

    public static final List<HmfTranscriptRegion> HLA_TRANSCRIPTS = Lists.newArrayList();

    public static final LociPosition LOCI_POSITION = new LociPosition();
}

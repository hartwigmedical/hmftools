package com.hartwig.hmftools.serve.vicc.filter;

import java.util.Set;

import com.google.common.collect.Sets;

final class FilterFactory {

    static final Set<String> FEATURE_KEYWORDS_TO_FILTER = Sets.newHashSet();

    static final Set<FilterKey> FEATURE_KEYS_TO_FILTER = Sets.newHashSet();

    static {
        populateFeatureNameKeywordsToFilter();
        populateFeatureKeysToFilter();
    }

    private static void populateFeatureNameKeywordsToFilter() {
        // We cannot determine methylation with WGS/WTS
        FEATURE_KEYWORDS_TO_FILTER.add("PROMOTER METHYLATION");
        FEATURE_KEYWORDS_TO_FILTER.add("PROMOTER DEMETHYLATION");
        FEATURE_KEYWORDS_TO_FILTER.add("Promoter Hypermethylation");
        FEATURE_KEYWORDS_TO_FILTER.add("METHYLATION");
        FEATURE_KEYWORDS_TO_FILTER.add("Hypermethylation");

        // We cannot determine epigenetic silencing with WGS/WTS
        FEATURE_KEYWORDS_TO_FILTER.add("Epigenetic Silencing");

        // We cannot determine phosphorylation with WGS/WTS
        FEATURE_KEYWORDS_TO_FILTER.add("PHOSPHORYLATION");

        // "Any" polymorphism is considered too vague.
        FEATURE_KEYWORDS_TO_FILTER.add("POLYMORPHISM");
        FEATURE_KEYWORDS_TO_FILTER.add("Single Nucleotide Polymorphism");

        // Copy number variation is too vague.
        FEATURE_KEYWORDS_TO_FILTER.add("COPY NUMBER VARIATION");
        FEATURE_KEYWORDS_TO_FILTER.add("EXPRESSION");
        FEATURE_KEYWORDS_TO_FILTER.add("expression");
        FEATURE_KEYWORDS_TO_FILTER.add("WILDTYPE");
        FEATURE_KEYWORDS_TO_FILTER.add("WILD TYPE");
        FEATURE_KEYWORDS_TO_FILTER.add("wild-type");
        FEATURE_KEYWORDS_TO_FILTER.add("wildtype");
        FEATURE_KEYWORDS_TO_FILTER.add("Wildtype");
        FEATURE_KEYWORDS_TO_FILTER.add("SERUM LEVELS");

        FEATURE_KEYWORDS_TO_FILTER.add("3' EXON DELETION");
        FEATURE_KEYWORDS_TO_FILTER.add("3' UTR MUTATION");
        FEATURE_KEYWORDS_TO_FILTER.add("3' UTR Polymorphism");
        FEATURE_KEYWORDS_TO_FILTER.add("5' TANDEM REPEAT");
        FEATURE_KEYWORDS_TO_FILTER.add("B2 DOMAIN MUTATION");
        FEATURE_KEYWORDS_TO_FILTER.add("CYTOPLASMIC MISLOCALIZATION");
        FEATURE_KEYWORDS_TO_FILTER.add("DNA binding domain deletions");
        FEATURE_KEYWORDS_TO_FILTER.add("DNA binding domain insertions");
        FEATURE_KEYWORDS_TO_FILTER.add("DNA binding domain missense mutations");
        FEATURE_KEYWORDS_TO_FILTER.add("DNA BINDING DOMAIN MUTATION");
        FEATURE_KEYWORDS_TO_FILTER.add("FLT3 internal tandem duplications");
        FEATURE_KEYWORDS_TO_FILTER.add("INTERNAL DUPLICATION");
        FEATURE_KEYWORDS_TO_FILTER.add("KINASE DOMAIN MUTATION");
        FEATURE_KEYWORDS_TO_FILTER.add("N-TERMINAL FRAME SHIFT");
        FEATURE_KEYWORDS_TO_FILTER.add("NUCLEAR TRANSLOCATION");
        FEATURE_KEYWORDS_TO_FILTER.add("PROMOTER MUTATION");
        FEATURE_KEYWORDS_TO_FILTER.add("Promoter Mutations");
        FEATURE_KEYWORDS_TO_FILTER.add("SH2 DOMAIN MUTATION");
        FEATURE_KEYWORDS_TO_FILTER.add("TKD MUTATION");
        FEATURE_KEYWORDS_TO_FILTER.add("EGFR CTD");
        FEATURE_KEYWORDS_TO_FILTER.add("DPYD*13 HOMOZYGOSITY");
        FEATURE_KEYWORDS_TO_FILTER.add("DPYD*2A HOMOZYGOSITY");
        FEATURE_KEYWORDS_TO_FILTER.add("S291fsX300");
        FEATURE_KEYWORDS_TO_FILTER.add("S70fsX93");
        FEATURE_KEYWORDS_TO_FILTER.add("COPY-NEUTRAL LOSS OF HETEROZYGOSITY");
        FEATURE_KEYWORDS_TO_FILTER.add("T599insTT");
        FEATURE_KEYWORDS_TO_FILTER.add("T76insTLDT");
        FEATURE_KEYWORDS_TO_FILTER.add("TA83del");
        FEATURE_KEYWORDS_TO_FILTER.add("V1790_A1996del");
        FEATURE_KEYWORDS_TO_FILTER.add("V555_L576del");
        FEATURE_KEYWORDS_TO_FILTER.add("BRAF V600E/K");
        FEATURE_KEYWORDS_TO_FILTER.add("D835H/Y");
        FEATURE_KEYWORDS_TO_FILTER.add("G12/G13");
        FEATURE_KEYWORDS_TO_FILTER.add("Q157P/R");
        FEATURE_KEYWORDS_TO_FILTER.add("S310F/Y");
        FEATURE_KEYWORDS_TO_FILTER.add("S893A/T");
        FEATURE_KEYWORDS_TO_FILTER.add("CASP8L");
        FEATURE_KEYWORDS_TO_FILTER.add("INTRON 6 MUTATION");
        FEATURE_KEYWORDS_TO_FILTER.add("SNP309");
        FEATURE_KEYWORDS_TO_FILTER.add("S34Y/F");

        FEATURE_KEYWORDS_TO_FILTER.add("RS11623866");
        FEATURE_KEYWORDS_TO_FILTER.add("RS16906252");
        FEATURE_KEYWORDS_TO_FILTER.add("RS2305035");
        FEATURE_KEYWORDS_TO_FILTER.add("RS2736100");
        FEATURE_KEYWORDS_TO_FILTER.add("RS3184504");
        FEATURE_KEYWORDS_TO_FILTER.add("RS34743033");
        FEATURE_KEYWORDS_TO_FILTER.add("RS34886328");
        FEATURE_KEYWORDS_TO_FILTER.add("RS3733542");
        FEATURE_KEYWORDS_TO_FILTER.add("RS3814960");
        FEATURE_KEYWORDS_TO_FILTER.add("RS4149056");
        FEATURE_KEYWORDS_TO_FILTER.add("RS461155");
        FEATURE_KEYWORDS_TO_FILTER.add("RS557806");
        FEATURE_KEYWORDS_TO_FILTER.add("RS61764370");
        FEATURE_KEYWORDS_TO_FILTER.add("RS67376798 HOMOZYGOSITY");
        FEATURE_KEYWORDS_TO_FILTER.add("rs681673");
        FEATURE_KEYWORDS_TO_FILTER.add("RS751402");

        // exon number not exist in the using transcript
        FEATURE_KEYWORDS_TO_FILTER.add("NPM1 EXON 12 MUTATION");
    }

    private static void populateFeatureKeysToFilter() {

    }

    private FilterFactory() {
    }
}

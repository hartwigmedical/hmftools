package com.hartwig.hmftools.serve.sources.vicc.filter;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

final class FilterFactory {

    static final Set<String> FEATURE_KEYWORDS_TO_FILTER = Sets.newHashSet();

    static final Set<FilterKey> FEATURE_KEYS_TO_FILTER = Sets.newHashSet();

    static {
        populateFeatureKeywordsToFilter();
        populateFeatureKeysToFilter();
    }

    private static void populateFeatureKeywordsToFilter() {
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


        // to specific
        FEATURE_KEYWORDS_TO_FILTER.add("EGFR (L858R,L861,G719,S768I)");
        FEATURE_KEYWORDS_TO_FILTER.add("ESR1 (E380Q,537,538,L536,P535H)");
        FEATURE_KEYWORDS_TO_FILTER.add("KIT (550-592,627-664,788-828,829-860)");
        FEATURE_KEYWORDS_TO_FILTER.add("MAP2K1 (Q56P,P124S,P124L;C121S)");
        FEATURE_KEYWORDS_TO_FILTER.add("PDGFRA (552-596,631-668,814-854)");
        FEATURE_KEYWORDS_TO_FILTER.add("FGFR2 (V565I,M536I,M538I,I548V,N550,E566G,L618M,K660E)");
        FEATURE_KEYWORDS_TO_FILTER.add("FLT3 (F691,D835,N676,Y842)");

        // exon number not exist in the using transcript //TODO move to vicc curation
        FEATURE_KEYWORDS_TO_FILTER.add("NPM1 EXON 12 MUTATION");

        // TODO All of below needs some further investigation & explanation
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
    }

    private static void populateFeatureKeysToFilter() {
        // Variants implying stop lost. They are real but not handled yet in SERVE (TODO DEV-1475)
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "VHL", "*214W (c.642A>G)"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "MLH1", "*757L"));

        // Synonymous variant, unlikely to have an impact
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "VHL", "R161R (c.481C>A)"));

        // Variant is fine, but interpretation currently not supported by SERVE
        // The unaligned delete is fine but the left-aligned insert lies outside of exon range (TODO DEV-1475)
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "KIT", "KIT K550_W557del "));

        // Variant is fine, but interpretation currently not supported by SERVE
        // The unaligned insert lies outside of exonic range, but the left-aligned insert is fine (TODO DEV-1475)
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "BRAF", "R506_K507insVLR"));

        // Variant is fine, but interpretation currently not supported by SERVE
        // The unaligned delete is fine but the left-aligned insert lies outside of exon range (TODO DEV-1475)
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "KIT", "K550_K558del"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "KIT", "K550_W557del"));
    }

    private FilterFactory() {
    }
}

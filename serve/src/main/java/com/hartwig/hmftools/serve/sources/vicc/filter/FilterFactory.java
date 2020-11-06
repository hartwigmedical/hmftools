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

        // We cannot determine the below events
        FEATURE_KEYWORDS_TO_FILTER.add("SERUM LEVELS");
        FEATURE_KEYWORDS_TO_FILTER.add("CYTOPLASMIC MISLOCALIZATION");
        FEATURE_KEYWORDS_TO_FILTER.add("NUCLEAR TRANSLOCATION");

        // "Any" polymorphism is considered too vague.
        FEATURE_KEYWORDS_TO_FILTER.add("POLYMORPHISM");
        FEATURE_KEYWORDS_TO_FILTER.add("Single Nucleotide Polymorphism");
        FEATURE_KEYWORDS_TO_FILTER.add("3' UTR Polymorphism");

        // Copy number variation is too vague.
        FEATURE_KEYWORDS_TO_FILTER.add("COPY NUMBER VARIATION");
        FEATURE_KEYWORDS_TO_FILTER.add("COPY-NEUTRAL LOSS OF HETEROZYGOSITY");

        // "Expression" is not observed on DNA level
        FEATURE_KEYWORDS_TO_FILTER.add("EXPRESSION");
        FEATURE_KEYWORDS_TO_FILTER.add("expression");

        // Wildtype evidence is ignored at this point.
        FEATURE_KEYWORDS_TO_FILTER.add("WILDTYPE");
        FEATURE_KEYWORDS_TO_FILTER.add("WILD TYPE");
        FEATURE_KEYWORDS_TO_FILTER.add("wild-type");
        FEATURE_KEYWORDS_TO_FILTER.add("wildtype");
        FEATURE_KEYWORDS_TO_FILTER.add("Wildtype");
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

        // exon 12 does not exist in canonical transcript of NPM1
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "NPM1", "EXON 12 MUTATION"));

        // We ignore variants specified in their RS identifier
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "FNTB", "RS11623866"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "MGMT", "RS16906252"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "CBLB", "RS2305035"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "TERT", "RS2736100"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "SH2B3", "RS3184504"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "TYMS", "RS34743033"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "MDM2", "RS34886328"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "KIT", "RS3733542"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "CDKN2A", "RS3814960"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "SLCO1B1", "RS4149056"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "ETS2", "RS461155"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "PPP1R15A", "RS557806"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "KRAS", "RS61764370"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "DPYD", "RS67376798 HOMOZYGOSITY"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "GADD45A", "rs681673"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "ERCC5", "RS751402"));
    }

    private FilterFactory() {
    }
}

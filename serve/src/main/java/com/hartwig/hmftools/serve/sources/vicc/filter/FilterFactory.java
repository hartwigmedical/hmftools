package com.hartwig.hmftools.serve.sources.vicc.filter;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

final class FilterFactory {

    static final Set<String> FEATURE_KEYWORDS_TO_FILTER = Sets.newHashSet();

    static final Set<String> FEATURES_TO_FILTER = Sets.newHashSet();

    static final Set<FilterKey> FEATURE_KEYS_TO_FILTER = Sets.newHashSet();

    static {
        populateFeatureKeywordsToFilter();
        populateFeaturesToFilter();
        populateFeatureKeysToFilter();
    }

    private static void populateFeatureKeywordsToFilter() {
        // Any wildtype evidence is ignored at this point.
        FEATURE_KEYWORDS_TO_FILTER.add("WILDTYPE");
        FEATURE_KEYWORDS_TO_FILTER.add("WILD TYPE");
        FEATURE_KEYWORDS_TO_FILTER.add("wild-type");
        FEATURE_KEYWORDS_TO_FILTER.add("wildtype");
        FEATURE_KEYWORDS_TO_FILTER.add("Wildtype");
    }

    private static void populateFeaturesToFilter() {
        // We cannot determine methylation with WGS/WTS
        FEATURES_TO_FILTER.add("METHYLATION");
        FEATURES_TO_FILTER.add("Hypermethylation");
        FEATURES_TO_FILTER.add("Promoter Hypermethylation");
        FEATURES_TO_FILTER.add("PROMOTER HYPERMETHYLATION");
        FEATURES_TO_FILTER.add("PROMOTER METHYLATION");
        FEATURES_TO_FILTER.add("PROMOTER DEMETHYLATION");

        // We cannot determine epigenetic silencing with WGS/WTS
        FEATURES_TO_FILTER.add("Epigenetic Silencing");

        // We cannot determine phosphorylation with WGS/WTS
        FEATURES_TO_FILTER.add("PHOSPHORYLATION");

        // We cannot determine the below events
        FEATURES_TO_FILTER.add("SERUM LEVELS");
        FEATURES_TO_FILTER.add("CYTOPLASMIC MISLOCALIZATION");
        FEATURES_TO_FILTER.add("NUCLEAR TRANSLOCATION");
        FEATURES_TO_FILTER.add("NUCLEAR EXPRESSION");

        // "Any" polymorphism is considered too vague.
        FEATURES_TO_FILTER.add("POLYMORPHISM");
        FEATURES_TO_FILTER.add("Single Nucleotide Polymorphism");
        FEATURES_TO_FILTER.add("3' UTR Polymorphism");

        // Copy number variation is too vague.
        FEATURES_TO_FILTER.add("COPY NUMBER VARIATION");
        FEATURES_TO_FILTER.add("COPY-NEUTRAL LOSS OF HETEROZYGOSITY");

        // "Expression" is not observed on DNA level
        FEATURES_TO_FILTER.add("CYTOPLASMIC EXPRESSION");
        FEATURES_TO_FILTER.add("ISOFORM EXPRESSION");
        FEATURES_TO_FILTER.add("EXPRESSION");
    }

    private static void populateFeatureKeysToFilter() {
        // Variants implying stop lost. They are real but not handled yet in SERVE
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "VHL", "*214W (c.642A>G)"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "MLH1", "*757L"));

        // Synonymous variant, unlikely to have an impact
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "VHL", "R161R (c.481C>A)"));

        // Variant is fine, but interpretation currently not supported by SERVE
        // The unaligned delete is fine but the left-aligned insert lies outside of exon range
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "KIT", "KIT K550_W557del "));

        // Variant is fine, but interpretation currently not supported by SERVE
        // The unaligned insert lies outside of exonic range, but the left-aligned insert is fine
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "BRAF", "R506_K507insVLR"));

        // Variant is fine, but interpretation currently not supported by SERVE
        // The unaligned delete is fine but the left-aligned insert lies outside of exon range
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

        // A similar type of filtered variant to RS variants.
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "MDM2", "SNP309"));

        // Event with this gene can be remove
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "BCL2", "IGH-BCL2"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "BCL2", "IGH-BCL2 Fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "CCND1", "IGH-CCND1 Fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "FGFR3", "IGH-FGFR3 Fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "MYC", "IGH-MYC Fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "NKX2-1", "IGH-NKX2 Fusion"));

        // We don't observe phosphorylation
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "EGFR", "Y1092 PHOSPHORYLATION"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "PRKAA2", "T172 PHOSPHORYLATION"));

        // We ignore generic "expression"
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CGI, "CYP17A1", "CYP17A1 expression"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "CDKN2A", "p16 EXPRESSION"));
    }

    private FilterFactory() {
    }
}

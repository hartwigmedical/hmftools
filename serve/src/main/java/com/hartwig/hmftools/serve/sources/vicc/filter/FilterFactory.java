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
        FEATURES_TO_FILTER.add("SH2 DOMAIN MUTATION");

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
        populateFilterForHotspotsThatCannotYetBeInterpreted();
        populateFilterForEventsInconsistentWithDriverCatalog();
        populateFilterForNonReportableFusions();

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

        // We don't observe phosphorylation
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "EGFR", "Y1092 PHOSPHORYLATION"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "PRKAA2", "T172 PHOSPHORYLATION"));

        // We ignore generic "expression"
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CGI, "CYP17A1", "CYP17A1 expression"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "CDKN2A", "p16 EXPRESSION"));
    }

    private static void populateFilterForHotspotsThatCannotYetBeInterpreted() {
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

        // Not sure if this fusion should be swapped, but if not then it would never be reported so can filter.
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CGI, "RET", "RET-TPCN1 fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CGI, "TPCN1", "RET-TPCN1 fusion"));
    }

    private static void populateFilterForNonReportableFusions() {
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CGI, "AKT3", "AKT3 fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CGI, "ERBB4", "ERBB4 fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CGI, "ESR1", "ESR1-YAP1 fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CGI, "NOTCH1", "NOTCH1 fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CGI, "NOTCH2", "NOTCH2 fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CGI, "YAP1", "ESR1-YAP1 fusion"));

        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "FOS", "TRUNCATING FUSION"));

        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "AKT2", "Fusions"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "AKT2", "BCAM-AKT2 Fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "AKT3", "Fusions"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "AKT3", "MAGI3-AKT3 Fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "ERBB4", "Fusions"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "ERBB4", "EZR-ERBB4 Fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "ESR1", "ESR1-CCDC170 Fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "ESR1", "ESR1-YAP1 Fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "ESR1", "Fusions"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "NKX2-1", "Fusions"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "NOTCH1", "Fusions"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "NOTCH1", "MIR143-NOTCH1 Fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "PIK3CB", "ACPP-PIK3CB Fusion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "PIK3CB", "Fusions"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "SEC16A", "SEC16A-NOTCH1 Fusion"));

        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "EGFR", "EGFR-PURB"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "NOTCH1", "NOTCH1  rearrange"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "PURB", "EGFR-PURB"));
    }

    private static void populateFilterForEventsInconsistentWithDriverCatalog() {
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CGI, "CDK12", "CDK12 amplification"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CGI, "EPHA2", "EPHA2 amplification"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CGI, "FLT1", "FLT1 overexpression"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CGI, "RB1", "RB1 overexpression"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CGI, "UGT1A1", "UGT1A1 biallelic inactivation"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CGI, "TPMT", "TPMT biallelic inactivation"));

        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "CCND3", "LOSS"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "FOXP1", "AMPLIFICATION"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "NOTCH1", "AMPLIFICATION"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "TYMS", "UNDEREXPRESSION"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "VEGFA", "UNDEREXPRESSION"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "JAK1", "OVEREXPRESSION"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.CIVIC, "RB1", "OVEREXPRESSION"));

        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "CBL", "CBL  del"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "CBL", "CBL  dec exp"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "BRAF", "BRAF  dec exp"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "ERBB2", "ERBB2  dec exp"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "EGFR", "EGFR  dec exp"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "FGFR1", "FGFR1  dec exp"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "FGFR2", "FGFR2  dec exp"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "FGFR3", "FGFR3  dec exp"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "MET", "MET  dec exp"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "PTPN11", "PTPN11  dec exp"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "SMO", "SMO  dec exp"));

        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "HRAS", "HRAS  inact mut"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "JAK2", "JAK2  inact mut"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "MET", "MET  negative"));

        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "CDKN2A", "CDKN2A  pos"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "CEBPA", "CEBPA  positive"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "NOTCH1", "NOTCH1  positive"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "PTEN", "PTEN  pos"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "PTEN", "PTEN  positive"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "RB1", "RB1  positive"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "TP53", "TP53  positive"));

        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "ATM", "ATM  over exp"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "CDH1", "CDH1  over exp"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "CDKN2A", "CDKN2A  over exp"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "GATA1", "GATA1  over exp"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "NOTCH1", "NOTCH1  act mut"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "NOTCH1", "NOTCH1  over exp"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.JAX, "VHL", "VHL  over exp"));

        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "KDM5C", "Overexpression"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "ERCC2", "Amplification"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "GATA3", "Amplification"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "RUNX1", "Amplification"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "NOTCH1", "Gain-of-function Mutations"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "FOXA1", "Deletion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "RAD21", "Deletion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "PTPN11", "Deletion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "CBL", "Deletion"));
        FEATURE_KEYS_TO_FILTER.add(new FilterKey(ViccSource.ONCOKB, "NPM1", "Deletion"));
    }

    private FilterFactory() {
    }
}

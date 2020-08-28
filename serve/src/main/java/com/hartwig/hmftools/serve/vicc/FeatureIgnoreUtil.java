package com.hartwig.hmftools.serve.vicc;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.jetbrains.annotations.NotNull;

public final class FeatureIgnoreUtil {

    private static final Set<String> FEATURE_KEYWORDS_TO_IGNORE = Sets.newHashSet();

    static {
        // We cannot determine methylation with WGS/WTS
        FEATURE_KEYWORDS_TO_IGNORE.add("PROMOTER METHYLATION");
        FEATURE_KEYWORDS_TO_IGNORE.add("PROMOTER DEMETHYLATION");
        FEATURE_KEYWORDS_TO_IGNORE.add("Promoter Hypermethylation");
        FEATURE_KEYWORDS_TO_IGNORE.add("METHYLATION");
        FEATURE_KEYWORDS_TO_IGNORE.add("Hypermethylation");

        // We cannot determine epigenetic silencing with WGS/WTS
        FEATURE_KEYWORDS_TO_IGNORE.add("Epigenetic Silencing");

        // We cannot determine phosphorylation with WGS/WTS
        FEATURE_KEYWORDS_TO_IGNORE.add("PHOSPHORYLATION");

        // "Any" polymorphism is considered too vague.
        FEATURE_KEYWORDS_TO_IGNORE.add("POLYMORPHISM");
        FEATURE_KEYWORDS_TO_IGNORE.add("Single Nucleotide Polymorphism");

        // Copy number variation is too vague.
        FEATURE_KEYWORDS_TO_IGNORE.add("COPY NUMBER VARIATION");
        FEATURE_KEYWORDS_TO_IGNORE.add("3' EXON DELETION"); //TODO determine to which event?
        FEATURE_KEYWORDS_TO_IGNORE.add("EXPRESSION");
        FEATURE_KEYWORDS_TO_IGNORE.add("expression");
        FEATURE_KEYWORDS_TO_IGNORE.add("WILDTYPE");
        FEATURE_KEYWORDS_TO_IGNORE.add("WILD TYPE");
        FEATURE_KEYWORDS_TO_IGNORE.add("wild-type");
        FEATURE_KEYWORDS_TO_IGNORE.add("wildtype");
        FEATURE_KEYWORDS_TO_IGNORE.add("Wildtype");
        FEATURE_KEYWORDS_TO_IGNORE.add("SERUM LEVELS");

        FEATURE_KEYWORDS_TO_IGNORE.add("3' UTR MUTATION");
        FEATURE_KEYWORDS_TO_IGNORE.add("3' UTR Polymorphism");
        FEATURE_KEYWORDS_TO_IGNORE.add("5' TANDEM REPEAT");
        FEATURE_KEYWORDS_TO_IGNORE.add("B2 DOMAIN MUTATION");
        FEATURE_KEYWORDS_TO_IGNORE.add("CYTOPLASMIC MISLOCALIZATION");
        FEATURE_KEYWORDS_TO_IGNORE.add("DNA binding domain deletions");
        FEATURE_KEYWORDS_TO_IGNORE.add("DNA binding domain insertions");
        FEATURE_KEYWORDS_TO_IGNORE.add("DNA binding domain missense mutations");
        FEATURE_KEYWORDS_TO_IGNORE.add("DNA BINDING DOMAIN MUTATION");
        FEATURE_KEYWORDS_TO_IGNORE.add("FLT3 internal tandem duplications");
        FEATURE_KEYWORDS_TO_IGNORE.add("INTERNAL DUPLICATION");
        FEATURE_KEYWORDS_TO_IGNORE.add("KINASE DOMAIN MUTATION");
        FEATURE_KEYWORDS_TO_IGNORE.add("N-TERMINAL FRAME SHIFT");
        FEATURE_KEYWORDS_TO_IGNORE.add("NUCLEAR TRANSLOCATION");
        FEATURE_KEYWORDS_TO_IGNORE.add("PROMOTER MUTATION");
        FEATURE_KEYWORDS_TO_IGNORE.add("Promoter Mutations");
        FEATURE_KEYWORDS_TO_IGNORE.add("SH2 DOMAIN MUTATION");
        FEATURE_KEYWORDS_TO_IGNORE.add("TKD MUTATION");
        FEATURE_KEYWORDS_TO_IGNORE.add("EGFR CTD");
        FEATURE_KEYWORDS_TO_IGNORE.add("DPYD*13 HOMOZYGOSITY");
        FEATURE_KEYWORDS_TO_IGNORE.add("DPYD*2A HOMOZYGOSITY");
        FEATURE_KEYWORDS_TO_IGNORE.add("S291fsX300");
        FEATURE_KEYWORDS_TO_IGNORE.add("S70fsX93");
        FEATURE_KEYWORDS_TO_IGNORE.add("COPY-NEUTRAL LOSS OF HETEROZYGOSITY");
        FEATURE_KEYWORDS_TO_IGNORE.add("T599insTT");
        FEATURE_KEYWORDS_TO_IGNORE.add("T76insTLDT");
        FEATURE_KEYWORDS_TO_IGNORE.add("TA83del");
        FEATURE_KEYWORDS_TO_IGNORE.add("V1790_A1996del");
        FEATURE_KEYWORDS_TO_IGNORE.add("V555_L576del");
        FEATURE_KEYWORDS_TO_IGNORE.add("BRAF V600E/K");
        FEATURE_KEYWORDS_TO_IGNORE.add("D835H/Y");
        FEATURE_KEYWORDS_TO_IGNORE.add("G12/G13");
        FEATURE_KEYWORDS_TO_IGNORE.add("Q157P/R");
        FEATURE_KEYWORDS_TO_IGNORE.add("S310F/Y");
        FEATURE_KEYWORDS_TO_IGNORE.add("S893A/T");
        FEATURE_KEYWORDS_TO_IGNORE.add("CASP8L");
        FEATURE_KEYWORDS_TO_IGNORE.add("INTRON 6 MUTATION");
        FEATURE_KEYWORDS_TO_IGNORE.add("SNP309");
        FEATURE_KEYWORDS_TO_IGNORE.add("S34Y/F");

        FEATURE_KEYWORDS_TO_IGNORE.add("RS11623866");
        FEATURE_KEYWORDS_TO_IGNORE.add("RS16906252");
        FEATURE_KEYWORDS_TO_IGNORE.add("RS2305035");
        FEATURE_KEYWORDS_TO_IGNORE.add("RS2736100");
        FEATURE_KEYWORDS_TO_IGNORE.add("RS3184504");
        FEATURE_KEYWORDS_TO_IGNORE.add("RS34743033");
        FEATURE_KEYWORDS_TO_IGNORE.add("RS34886328");
        FEATURE_KEYWORDS_TO_IGNORE.add("RS3733542");
        FEATURE_KEYWORDS_TO_IGNORE.add("RS3814960");
        FEATURE_KEYWORDS_TO_IGNORE.add("RS4149056");
        FEATURE_KEYWORDS_TO_IGNORE.add("RS461155");
        FEATURE_KEYWORDS_TO_IGNORE.add("RS557806");
        FEATURE_KEYWORDS_TO_IGNORE.add("RS61764370");
        FEATURE_KEYWORDS_TO_IGNORE.add("RS67376798 HOMOZYGOSITY");
        FEATURE_KEYWORDS_TO_IGNORE.add("rs681673");
        FEATURE_KEYWORDS_TO_IGNORE.add("RS751402");

        // exon number not exist in the using transcript
        FEATURE_KEYWORDS_TO_IGNORE.add("NPM1 EXON 12 MUTATION");

    }

    private FeatureIgnoreUtil() {
    }

    public static boolean canIgnore(@NotNull Feature feature) {
        for (String ignoreKeyword : FEATURE_KEYWORDS_TO_IGNORE) {
            if (feature.name().contains(ignoreKeyword)) {
                return true;
            } else if (feature.description() != null){
                if (feature.description().contains(ignoreKeyword)) {
                    return true;
                }
            }
        }

        return false;
    }
}

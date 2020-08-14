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
        FEATURE_KEYWORDS_TO_IGNORE.add("WILDTYPE");
        FEATURE_KEYWORDS_TO_IGNORE.add("wild-type");
        FEATURE_KEYWORDS_TO_IGNORE.add("wildtype");

    }

    private FeatureIgnoreUtil() {
    }

    public static boolean canIgnore(@NotNull Feature feature) {
        for (String ignoreKeyword : FEATURE_KEYWORDS_TO_IGNORE) {
            if (feature.name().contains(ignoreKeyword)) {
                return true;
            }
        }

        return false;
    }
}

package com.hartwig.hmftools.vicc.annotation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.ExclusiveEventClassifier;

import org.jetbrains.annotations.NotNull;

public class GeneLevelClassifier implements EventClassifier {

    public static final Set<String> GENERIC_GENE_LEVEL_KEYWORDS = Sets.newHashSet("MUTATION",
            "mutant",
            "mut",
            "TRUNCATING MUTATION",
            "Truncating Mutations",
            "feature_truncation",
            "FRAMESHIFT TRUNCATION",
            "FRAMESHIFT MUTATION",
            "ALTERATION");

    public static final Set<String> INACTIVATING_GENE_LEVEL_KEYWORDS = Sets.newHashSet("inact mut",
            "biallelic inactivation",
            "Loss Of Function Variant",
            "Loss Of Heterozygosity",
            "DELETERIOUS MUTATION",
            "negative",
            "BIALLELIC INACTIVATION",
            "LOSS-OF-FUNCTION");

    public static final Set<String> ACTIVATING_GENE_LEVEL_KEYWORDS = Sets.newHashSet("Gain-of-function Mutations",
            "Gain-of-Function",
            "act mut",
            "ACTIVATING MUTATION",
            "Oncogenic Mutations",
            "pos",
            "positive",
            "oncogenic mutation");

    private static final String EXON_KEYWORD = "exon";

    @NotNull
    public static EventClassifier create(@NotNull List<EventClassifier> excludingEventClassifiers) {
        return new ExclusiveEventClassifier(excludingEventClassifiers, new GeneLevelClassifier());
    }

    private GeneLevelClassifier() {
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        if (event.toLowerCase().contains(EXON_KEYWORD)) {
            return false;
        }

        for (String keyword : GENERIC_GENE_LEVEL_KEYWORDS) {
            if (event.contains(keyword)) {
                return true;
            }
        }

        for (String keyword : INACTIVATING_GENE_LEVEL_KEYWORDS) {
            if (event.contains(keyword)) {
                return true;
            }
        }

        for (String keyword : ACTIVATING_GENE_LEVEL_KEYWORDS) {
            if (event.contains(keyword)) {
                return true;
            }
        }

        return event.trim().equals(gene);
    }
}

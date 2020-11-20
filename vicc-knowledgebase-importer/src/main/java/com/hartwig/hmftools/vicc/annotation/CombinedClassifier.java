package com.hartwig.hmftools.vicc.annotation;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventMatcher;

import org.jetbrains.annotations.NotNull;

class CombinedClassifier implements EventMatcher {

    private static final Map<String, Set<String>> COMBINED_EVENTS_PER_GENE = Maps.newHashMap();

    static {
        COMBINED_EVENTS_PER_GENE.put("EGFR", Sets.newHashSet("Ex19 del L858R"));
        COMBINED_EVENTS_PER_GENE.put("BRAF", Sets.newHashSet("p61BRAF-V600E", "V600E AMPLIFICATION"));
    }

    public CombinedClassifier() {
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        Set<String> entriesForGene = COMBINED_EVENTS_PER_GENE.get(gene);
        if (entriesForGene != null) {
            if (entriesForGene.contains(event)) {
                return true;
            }
        }

        if ((event.contains(",") || event.contains(";")) && !event.toLowerCase().contains(" or ")) {
            return true;
        } else if (event.contains("+") && !event.toLowerCase().contains("c.") && !event.contains(">")) {
            return true;
        } else if (event.contains("/")) {
            return false;
        } else if (event.trim().contains(" ")) {
            String[] parts = event.trim().replace("  ", " ").split(" ");
            if (parts[0].contains("-")) {
                // Hotspots or amplifications on fusion genes are considered combined.
                return HotspotClassifier.isValidProteinAnnotation(parts[1]) || AmplificationClassifier.isTypicalAmplification(parts[1]);
            }
        }

        return false;
    }
}

package com.hartwig.hmftools.vicc.annotation;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class CombinedClassifier {

    private static final Map<String, Set<String>> FUSION_PAIR_AND_EXON_RANGES_PER_GENE = Maps.newHashMap();

    private static final Map<String, Set<String>> COMBINED_EVENTS_PER_GENE = Maps.newHashMap();

    static {
        Set<String> kitSet = Sets.newHashSet("EXON 11 MUTATION", "Exon 11 mutations", "Exon 11 deletions");
        Set<String> metSet = Sets.newHashSet("EXON 14 SKIPPING MUTATION");

        FUSION_PAIR_AND_EXON_RANGES_PER_GENE.put("KIT", kitSet);
        FUSION_PAIR_AND_EXON_RANGES_PER_GENE.put("MET", metSet);

        COMBINED_EVENTS_PER_GENE.put("EGFR", Sets.newHashSet("Ex19 del L858R"));
        COMBINED_EVENTS_PER_GENE.put("BRAF", Sets.newHashSet("p61BRAF-V600E"));
    }

    private CombinedClassifier() {
    }

    public static boolean isFusionPairAndGeneRangeExon(@NotNull String featureName, @Nullable String gene) {
        Set<String> entries = FUSION_PAIR_AND_EXON_RANGES_PER_GENE.get(gene);
        if (entries != null) {
            return entries.contains(featureName);
        }

        return false;
    }

    public static boolean isCombinedEvent(@NotNull String featureName, @Nullable String gene) {
        Set<String> entriesForGene = COMBINED_EVENTS_PER_GENE.get(gene);
        if (entriesForGene != null) {
            if (entriesForGene.contains(featureName)) {
                return true;
            }
        }

        if ((featureName.contains(",") || featureName.contains(";")) && !featureName.toLowerCase().contains(" or ")) {
            return true;
        } else if (featureName.contains("+") && !featureName.toLowerCase().contains("c.") && !featureName.contains(">")) {
            return true;
        } else if (featureName.contains("/")) {
            return false;
        } else if (featureName.trim().contains(" ")) {
            String[] parts = featureName.trim().replace("  ", " ").split(" ");
            if (parts[0].contains("-")) {
                // Hotspots or amplifications on fusion genes are considered combined.
                return HotspotClassifier.isHotspot(parts[1]) || CopyNumberClassifier.isAmplification(parts[1], gene);
            }
        }

        return false;
    }
}

package com.hartwig.hmftools.vicc.annotation;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

final class CombinedClassifier {

    private static final Map<String, Set<String>> FUSION_PAIR_AND_EXON_RANGES_PER_GENE = Maps.newHashMap();

    static {
        //TODO: check if more EXON_DEL_DUP fusions need to be added
        Set<String> kitSet = Sets.newHashSet("EXON 11 MUTATION", "Exon 11 mutations", "Exon 11 deletions");
        Set<String> metSet = Sets.newHashSet("EXON 14 SKIPPING MUTATION");

        FUSION_PAIR_AND_EXON_RANGES_PER_GENE.put("KIT", kitSet);
        FUSION_PAIR_AND_EXON_RANGES_PER_GENE.put("MET", metSet);
    }

    private CombinedClassifier() {
    }

    public static boolean isFusionPairAndGeneRangeExon(@NotNull String featureName, @NotNull String gene) {
        Set<String> entries = FUSION_PAIR_AND_EXON_RANGES_PER_GENE.get(gene);
        if (entries != null) {
            return entries.contains(featureName);
        }

        return false;
    }

    public static boolean isCombinedEvent(@NotNull String featureName) {
        if (featureName.contains("+") && !featureName.toLowerCase().contains("c.") && !featureName.contains(">")) {
            return true;
        } else if (featureName.contains("insertion")) {
            int countInsertion = featureName.split("insertion").length - 1;
            return countInsertion > 1;
        } else if (featureName.contains("deletion")) {
            int countDeletion = featureName.split("deletion").length - 1;
            return countDeletion > 1;
        } else if (featureName.contains("frameshift")) {
            int countFrameshift = featureName.split("frameshift").length - 1;
            return countFrameshift > 1;
        } else if (featureName.contains("insertions") && featureName.contains("deletion")) {
            int countCombined = (featureName.split("insertion").length - 1) + (featureName.split("deletion").length - 1);
            return countCombined > 1;
        } else if (featureName.contains("splice")) {
            int countSplice = featureName.split("splice").length - 1;
            return countSplice > 1;
        }

        return featureName.equals("p61BRAF-V600E");
    }
}

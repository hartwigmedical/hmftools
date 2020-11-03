package com.hartwig.hmftools.vicc.annotation;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class CombinedClassifier {

    //TODO: check if more EXON_DEL_DUP fusions need to be added
    private static final Set<String> FUSION_PAIR_AND_EXON_RANGE =
            Sets.newHashSet("KIT EXON 11 MUTATION", "KIT Exon 11 mutations", "KIT Exon 11 deletions", "MET EXON 14 SKIPPING MUTATION");

    private CombinedClassifier() {
    }

    public static boolean isFusionPairAndGeneRangeExon(@Nullable String featureDescription) {
        return FUSION_PAIR_AND_EXON_RANGE.contains(featureDescription);
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

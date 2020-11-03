package com.hartwig.hmftools.vicc.annotation;

import org.jetbrains.annotations.NotNull;

public final class CombinedClassifier {

    private CombinedClassifier() {
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

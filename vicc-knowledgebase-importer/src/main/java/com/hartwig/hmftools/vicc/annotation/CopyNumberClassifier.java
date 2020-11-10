package com.hartwig.hmftools.vicc.annotation;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class CopyNumberClassifier {

    private static final Set<String> AMPLIFICATION_KEYWORDS =
            Sets.newHashSet("Amplification", "amplification", "AMPLIFICATION", "amp", "overexpression", "OVEREXPRESSION", "Overexpression");

    private static final Set<String> AMPLIFICATION_KEY_PHRASES = Sets.newHashSet("over exp");

    private static final Set<String> DELETION_KEYWORDS =
            Sets.newHashSet("Deletion", "deletion", "DELETION", "del", "undexpression", "UNDEREXPRESSION", "loss", "LOSS");

    private static final Set<String> DELETION_KEY_PHRASES = Sets.newHashSet("dec exp", "Copy Number Loss");

    private static final Set<String> KEYWORDS_TO_SKIP_FOR_DELETION = Sets.newHashSet("exon", "EXON", "Exon", "Ex19", "inframe");

    private CopyNumberClassifier() {
    }

    public static boolean isAmplification(@NotNull String featureName, @Nullable String gene) {
        if (CombinedClassifier.isCombinedEvent(featureName, gene)) {
            return false;
        }

        return extractAmplificationDeletionType(featureName) == CopyNumberType.AMPLIFICATION;
    }

    public static boolean isDeletion(@NotNull String featureName, @Nullable String gene) {
        if (CombinedClassifier.isCombinedEvent(featureName, gene)) {
            return false;
        } else {
            for (String skipTerm : KEYWORDS_TO_SKIP_FOR_DELETION) {
                if (featureName.contains(skipTerm)) {
                    return false;
                }
            }

            return extractAmplificationDeletionType(featureName) == CopyNumberType.DELETION;
        }
    }

    @Nullable
    private static CopyNumberType extractAmplificationDeletionType(@NotNull String featureName) {
        String[] words = featureName.split(" ");
        for (String keyword : AMPLIFICATION_KEYWORDS) {
            for (String word : words) {
                if (word.equals(keyword)) {
                    return CopyNumberType.AMPLIFICATION;
                }
            }
        }

        for (String keyword : DELETION_KEYWORDS) {
            for (String word : words) {
                if (word.equals(keyword)) {
                    return CopyNumberType.DELETION;
                }
            }
        }

        for (String keyPhrase : AMPLIFICATION_KEY_PHRASES) {
            if (featureName.contains(keyPhrase)) {
                return CopyNumberType.AMPLIFICATION;
            }
        }

        for (String keyPhrase : DELETION_KEY_PHRASES) {
            if (featureName.contains(keyPhrase)) {
                return CopyNumberType.DELETION;
            }
        }

        return null;
    }

    private enum CopyNumberType {
        AMPLIFICATION,
        DELETION
    }
}

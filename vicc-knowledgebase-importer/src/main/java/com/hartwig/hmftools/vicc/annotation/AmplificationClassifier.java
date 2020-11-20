package com.hartwig.hmftools.vicc.annotation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.ExclusiveEventClassifier;

import org.jetbrains.annotations.NotNull;

class AmplificationClassifier implements EventClassifier {

    private static final Set<String> AMPLIFICATION_KEYWORDS =
            Sets.newHashSet("Amplification", "amplification", "AMPLIFICATION", "amp", "overexpression", "OVEREXPRESSION", "Overexpression");

    private static final Set<String> AMPLIFICATION_KEY_PHRASES = Sets.newHashSet("over exp");

    @NotNull
    public static EventClassifier create(@NotNull List<EventClassifier> excludingEventClassifiers) {
        return new ExclusiveEventClassifier(excludingEventClassifiers, new AmplificationClassifier());
    }

    private AmplificationClassifier() {
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        return isTypicalAmplification(event);
    }

    public static boolean isTypicalAmplification(@NotNull String event) {
        String[] words = event.split(" ");
        for (String keyword : AMPLIFICATION_KEYWORDS) {
            for (String word : words) {
                if (word.equals(keyword)) {
                    return true;
                }
            }
        }

        for (String keyPhrase : AMPLIFICATION_KEY_PHRASES) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }

        return false;
    }
}

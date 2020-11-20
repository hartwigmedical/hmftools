package com.hartwig.hmftools.vicc.annotation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventMatcher;
import com.hartwig.hmftools.common.serve.classification.ExclusiveEventMatcher;

import org.jetbrains.annotations.NotNull;

class AmplificationClassifier implements EventMatcher {

    private static final Set<String> AMPLIFICATION_KEYWORDS =
            Sets.newHashSet("Amplification", "amplification", "AMPLIFICATION", "amp", "overexpression", "OVEREXPRESSION", "Overexpression");

    private static final Set<String> AMPLIFICATION_KEY_PHRASES = Sets.newHashSet("over exp");

    @NotNull
    public static EventMatcher create(@NotNull List<EventMatcher> excludingEventMatchers) {
        return new ExclusiveEventMatcher(excludingEventMatchers, new AmplificationClassifier());
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

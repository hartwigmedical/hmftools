package com.hartwig.hmftools.vicc.annotation;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ComplexClassifier {

    private static final Map<String, Set<String>> COMPLEX_EVENTS_PER_GENE = Maps.newHashMap();

    static {
        Set<String> alkSet = Sets.newHashSet("ALK inframe insertion (1151T)");
        COMPLEX_EVENTS_PER_GENE.put("ALK", alkSet);
    }

    public static boolean isComplexEvent(@NotNull String featureName, @Nullable String gene) {
        Set<String> entriesForGene = COMPLEX_EVENTS_PER_GENE.get(gene);
        if (entriesForGene != null) {
            return !entriesForGene.contains(featureName);
        }
        return false;
    }
}

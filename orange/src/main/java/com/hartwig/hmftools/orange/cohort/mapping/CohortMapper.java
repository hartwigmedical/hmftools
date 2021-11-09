package com.hartwig.hmftools.orange.cohort.mapping;

import java.util.List;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.doid.DoidParents;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class CohortMapper {

    @NotNull
    private final DoidParents doidParentModel;
    @NotNull
    private final List<CohortMapping> mappings;

    public CohortMapper(@NotNull final DoidParents doidParentModel, @NotNull final List<CohortMapping> mappings) {
        this.doidParentModel = doidParentModel;
        this.mappings = mappings;
    }

    @NotNull
    public String cancerTypeForDoids(@NotNull Set<String> doids) {
        Multimap<String, CohortMapping> positiveMatchesPerDoid = ArrayListMultimap.create();

        for (String doid : doids) {
            Set<String> expandedDoids = Sets.newHashSet();

            expandedDoids.add(doid);
            expandedDoids.addAll(doidParentModel.parents(doid));

            for (CohortMapping mapping : mappings) {
                if (isMatch(mapping, expandedDoids)) {
                    positiveMatchesPerDoid.put(doid, mapping);
                }
            }
        }

        return Strings.EMPTY;
    }

    private static boolean isMatch(@NotNull CohortMapping mapping, @NotNull Set<String> expandedDoids) {
        boolean include = false;
        for (String doid : mapping.include()) {
            if (expandedDoids.contains(doid)) {
                include = true;
                break;
            }
        }

        boolean exclude = false;
        for (String doid : mapping.exclude()) {
            if (expandedDoids.contains(doid)) {
                exclude = true;
                break;
            }
        }

        return include && !exclude;
    }
}

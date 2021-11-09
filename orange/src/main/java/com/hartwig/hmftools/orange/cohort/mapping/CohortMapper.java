package com.hartwig.hmftools.orange.cohort.mapping;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.doid.DoidParents;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CohortMapper {

    private static final Logger LOGGER = LogManager.getLogger(CohortMapper.class);

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
            for (CohortMapping mapping : mappings) {
                if (isMatch(mapping, doid, doidParentModel.parents(doid))) {
                    positiveMatchesPerDoid.put(doid, mapping);
                }
            }
        }

        if (positiveMatchesPerDoid.isEmpty()) {
            LOGGER.warn("No positive doid matches found for doids '{}'. Reverting to cohort '{}'", doids, CohortConstants.COHORT_OTHER);
            return CohortConstants.COHORT_OTHER;
        } else {
            return pickBestCancerType(positiveMatchesPerDoid);
        }
    }

    private static boolean isMatch(@NotNull CohortMapping mapping, @NotNull String child, @NotNull Set<String> parents) {
        boolean include = false;
        for (String doid : mapping.include()) {
            if (child.equals(doid)) {
                include = true;
                break;
            } else if (parents.contains(doid) && mapping.rule() != MappingRule.EXACT_MATCH) {
                include = true;
                break;
            }
        }

        boolean exclude = false;
        for (String doid : mapping.exclude()) {
            if (parents.contains(doid) || child.equals(doid)) {
                exclude = true;
                break;
            }
        }

        return include && !exclude;
    }

    @NotNull
    private static String pickBestCancerType(@NotNull Multimap<String, CohortMapping> positiveMatchesPerDoid) {
        List<CohortMapping> bestMappings = Lists.newArrayList();
        for (Map.Entry<String, Collection<CohortMapping>> entry : positiveMatchesPerDoid.asMap().entrySet()) {
            Collection<CohortMapping> mappings = entry.getValue();
            if (mappings.size() == 1) {
                bestMappings.add(mappings.iterator().next());
            } else if (mappings.size() > 1) {
                LOGGER.warn("DOID '{}' led to multiple mappings: {}. Reverting to cohort '{}'",
                        entry.getKey(),
                        mappings,
                        CohortConstants.COHORT_OTHER);
                return CohortConstants.COHORT_OTHER;
            }
        }

        bestMappings.sort(new PreferenceRankComparator());
        return bestMappings.get(0).cancerType();
    }

    private static class PreferenceRankComparator implements Comparator<CohortMapping> {

        @Override
        public int compare(@NotNull CohortMapping mapping1, @NotNull CohortMapping mapping2) {
            if (mapping1.preferenceRank() == mapping2.preferenceRank()) {
                return 0;
            } else {
                return mapping1.preferenceRank() > mapping2.preferenceRank() ? 1 : -1;
            }
        }
    }

}

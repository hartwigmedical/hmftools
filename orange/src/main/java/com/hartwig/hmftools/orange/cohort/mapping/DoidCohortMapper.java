package com.hartwig.hmftools.orange.cohort.mapping;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.orange.cohort.datamodel.Sample;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class DoidCohortMapper implements CohortMapper {

    private static final Logger LOGGER = LogManager.getLogger(DoidCohortMapper.class);

    @NotNull
    private final DoidParents doidParentModel;
    @NotNull
    private final List<CohortMapping> mappings;

    public DoidCohortMapper(@NotNull final DoidParents doidParentModel, @NotNull final List<CohortMapping> mappings) {
        this.doidParentModel = doidParentModel;
        this.mappings = mappings;
    }

    @Nullable
    @Override
    public String cancerTypeForSample(@NotNull Sample sample) {
        Multimap<String, CohortMapping> positiveMatchesPerDoid = ArrayListMultimap.create();

        if (CohortConstants.DOID_COMBINATIONS_TO_MAP_TO_OTHER.contains(sample.doids())) {
            LOGGER.debug("Mapping {} to {} because of specific doid combination: {}",
                    sample.sampleId(),
                    CohortConstants.COHORT_OTHER,
                    sample.doids());
            return CohortConstants.COHORT_OTHER;
        }

        for (String doid : sample.doids()) {
            for (CohortMapping mapping : mappings) {
                if (isMatch(mapping, doid, doidParentModel.parents(doid))) {
                    positiveMatchesPerDoid.put(doid, mapping);
                }
            }
        }

        if (positiveMatchesPerDoid.isEmpty()) {
            LOGGER.warn("No DOID matches found for {}", toString(sample));
            return null;
        } else {
            return pickBestCancerType(sample, positiveMatchesPerDoid);
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

    @Nullable
    private static String pickBestCancerType(@NotNull Sample sample, @NotNull Multimap<String, CohortMapping> positiveMatchesPerDoid) {
        List<CohortMapping> bestMappings = Lists.newArrayList();
        for (Map.Entry<String, Collection<CohortMapping>> entry : positiveMatchesPerDoid.asMap().entrySet()) {
            Collection<CohortMapping> mappings = entry.getValue();
            if (mappings.size() == 1) {
                bestMappings.add(mappings.iterator().next());
            } else if (mappings.size() > 1) {
                String doid = entry.getKey();
                LOGGER.warn("DOID '{}' for {} matched to multiple cancer types: '{}'", doid, sample.sampleId(), toString(mappings));
                return null;
            }
        }

        bestMappings.sort(new PreferenceRankComparator());
        if (bestMappings.size() > 1) {
            CohortMapping bestMapping = bestMappings.get(0);
            for (int i = 1; i < bestMappings.size(); i++) {
                CohortMapping compare = bestMappings.get(i);
                if (bestMapping.preferenceRank() == compare.preferenceRank() && !bestMapping.cancerType().equals(compare.cancerType())) {
                    LOGGER.warn("Multiple different cancer types for {} with same preference rank: '{}'",
                            toString(sample),
                            toString(bestMappings));
                    return null;
                }
            }
        }

        return bestMappings.get(0).cancerType();
    }

    @NotNull
    private static String toString(@NotNull Collection<CohortMapping> mappings) {
        StringJoiner joiner = new StringJoiner(", ");
        for (CohortMapping mapping : mappings) {
            joiner.add(mapping.cancerType());
        }

        return joiner.toString();
    }

    @NotNull
    private static String toString(@NotNull Sample sample) {
        StringJoiner joiner = new StringJoiner(", ");
        for (String doid : sample.doids()) {
            joiner.add(doid);
        }
        return sample.sampleId() + " (doids=" + joiner + ")";
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

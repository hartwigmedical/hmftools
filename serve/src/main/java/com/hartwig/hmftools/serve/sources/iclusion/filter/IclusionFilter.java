package com.hartwig.hmftools.serve.sources.iclusion.filter;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.iclusion.data.IclusionMutation;
import com.hartwig.hmftools.iclusion.data.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.data.IclusionMutationLogicType;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionMutationCondition;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionTrial;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class IclusionFilter {

    private static final Logger LOGGER = LogManager.getLogger(IclusionFilter.class);

    @NotNull
    private final Set<String> filteredMutations = Sets.newHashSet();
    @NotNull
    private final Set<FilterKey> filteredKeys = Sets.newHashSet();

    public IclusionFilter() {
    }

    @NotNull
    public List<IclusionTrial> run(@NotNull List<IclusionTrial> trials) {
        // At this point we filter the following types of trials and conditions:
        //  - Trials without molecular inclusion criteria
        //  - Trials with empty acronyms.
        //  - Conditions that are grouped with a logic type other than OR (no combined events supported yet)
        //  - Conditions that are negated (not yet supported).

        List<IclusionTrial> filteredTrials = Lists.newArrayList();
        for (IclusionTrial trial : trials) {
            if (!trial.acronym().isEmpty()) {
                List<IclusionMutationCondition> filteredConditions = Lists.newArrayList();

                for (IclusionMutationCondition condition : trial.mutationConditions()) {
                    List<IclusionMutation> filteredMutations = Lists.newArrayList();
                    if (condition.logicType() == IclusionMutationLogicType.OR) {
                        for (IclusionMutation mutation : condition.mutations()) {
                            if (!mutation.negation()) {
                                if (include(mutation)) {
                                    filteredMutations.add(mutation);
                                } else {
                                    LOGGER.debug("Filtering mutation from {}: '{}' on '{}'",
                                            trial.acronym(),
                                            mutation.name(),
                                            mutation.gene());
                                }
                            } else {
                                LOGGER.debug("Filtering negated mutation from {}: '{}'", trial.acronym(), mutation);
                            }
                        }

                        if (!filteredMutations.isEmpty()) {
                            filteredConditions.add(ImmutableIclusionMutationCondition.builder()
                                    .from(condition)
                                    .mutations(filteredMutations)
                                    .build());
                        }
                    } else {
                        LOGGER.debug("Filtering non-OR mutation condition from {}: '{}'", trial.acronym(), condition);
                    }
                }

                if (!filteredConditions.isEmpty()) {
                    filteredTrials.add(ImmutableIclusionTrial.builder().from(trial).mutationConditions(filteredConditions).build());
                }
            } else {
                LOGGER.debug("Filtering trial for missing an acronym: '{}'", trial);
            }
        }

        return filteredTrials;
    }

    public void reportUnusedFilterEntries() {
        int unusedMutationCount = 0;
        for (String mutation : FilterFactory.MUTATIONS_TO_FILTER) {
            if (!filteredMutations.contains(mutation)) {
                unusedMutationCount++;
                LOGGER.warn("Mutation '{}' hasn't been used for iClusion filtering", mutation);
            }
        }

        int unusedKeyCount = 0;
        for (FilterKey key : FilterFactory.MUTATION_KEYS_TO_FILTER) {
            if (!filteredKeys.contains(key)) {
                unusedKeyCount++;
                LOGGER.warn("Key '{}' hasn't been used for iClusion filtering", key);
            }
        }

        LOGGER.debug("Found {} unused iClusion filter mutations and {} unused iClusion filter keys",
                unusedMutationCount,
                unusedKeyCount);
    }

    @VisibleForTesting
    boolean include(@NotNull IclusionMutation mutation) {
        for (String mutationToFilter : FilterFactory.MUTATIONS_TO_FILTER) {
            if (mutation.name().equals(mutationToFilter)) {
                filteredMutations.add(mutationToFilter);
                return false;
            }
        }

        FilterKey filterKey = new FilterKey(mutation.gene(), mutation.name());
        if (FilterFactory.MUTATION_KEYS_TO_FILTER.contains(filterKey)) {
            filteredKeys.add(filterKey);
            return false;
        }

        return true;
    }
}

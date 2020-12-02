package com.hartwig.hmftools.serve.sources.iclusion;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.iclusion.data.IclusionMutation;
import com.hartwig.hmftools.iclusion.data.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.serve.sources.ExtractionOutput;
import com.hartwig.hmftools.serve.sources.ImmutableExtractionOutput;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class IclusionExtractor {

    private static final Logger LOGGER = LogManager.getLogger(IclusionExtractor.class);

    @NotNull
    public ExtractionOutput extractFromIclusionTrials(@NotNull List<IclusionTrial> trials) {
        // We assume filtered trials (no empty acronyms, only OR mutations, and no negated mutations

        Set<String> uniqueCancerTypes = Sets.newTreeSet();
        for (IclusionTrial trial : trials) {
            List<ActionableTrial> actionableTrials = ActionableTrialFactory.toActionableTrials(trial);
            for (ActionableTrial actionableTrial : actionableTrials) {
                LOGGER.debug("Generated {} based off {}", actionableTrial, trial);
                uniqueCancerTypes.add(actionableTrial.cancerType());
            }

            for (IclusionMutationCondition mutationCondition : trial.mutationConditions()) {
                for (IclusionMutation mutation : mutationCondition.mutations()) {
                    LOGGER.debug("Interpreting '{}' on '{}' for {}", mutation.name(), mutation.gene(), trial.acronym());
                }
            }
        }

        LOGGER.debug("Printing {} unique cancer types", uniqueCancerTypes.size());
        for (String cancerType : uniqueCancerTypes) {
            LOGGER.debug(" - {}", cancerType);
        }

        return ImmutableExtractionOutput.builder().build();
    }
}

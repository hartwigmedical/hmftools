package com.hartwig.hmftools.serve.sources.iclusion;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutation;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.serve.ExtractionResult;
import com.hartwig.hmftools.serve.ImmutableExtractionResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class IclusionExtractor {

    private static final Logger LOGGER = LogManager.getLogger(IclusionExtractor.class);

    @NotNull
    public ExtractionResult extractFromIclusionTrials(@NotNull List<IclusionTrial> trials) {
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

        // TODO Implement!
        return ImmutableExtractionResult.builder().build();
    }
}

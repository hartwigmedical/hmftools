package com.hartwig.hmftools.serve.sources.iclusion;

import java.util.List;

import com.hartwig.hmftools.iclusion.datamodel.IclusionMutation;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.serve.extraction.EventExtractor;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ImmutableExtractionResult;
import com.hartwig.hmftools.serve.util.ProgressTracker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class IclusionExtractor {

    private static final Logger LOGGER = LogManager.getLogger(IclusionExtractor.class);

    @NotNull
    private final EventExtractor eventExtractor;

    public IclusionExtractor(@NotNull final EventExtractor eventExtractor) {
        this.eventExtractor = eventExtractor;
    }

    @NotNull
    public ExtractionResult extract(@NotNull List<IclusionTrial> trials) {
        // We assume filtered trials (no empty acronyms, only OR mutations, and no negated mutations

        ProgressTracker tracker = new ProgressTracker("iClusion", trials.size());
        for (IclusionTrial trial : trials) {
            List<ActionableTrial> actionableTrials = ActionableTrialFactory.toActionableTrials(trial);
            for (ActionableTrial actionableTrial : actionableTrials) {
                LOGGER.debug("Generated {} based off {}", actionableTrial, trial);
            }

            for (IclusionMutationCondition mutationCondition : trial.mutationConditions()) {
                for (IclusionMutation mutation : mutationCondition.mutations()) {
                    LOGGER.debug("Interpreting '{}' on '{}' for {}", mutation.name(), mutation.gene(), trial.acronym());
                }
            }

            tracker.update();
        }


        // TODO Implement!
        return ImmutableExtractionResult.builder().build();
    }
}
